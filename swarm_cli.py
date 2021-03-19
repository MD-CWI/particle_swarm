#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import tempfile
import shutil
import os
import sys
import math
import numpy as np
import pandas as pd
from subprocess import check_output


def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Command line interface to compute electron
        transport data from swarm simulations.
        Author: Jannis Teunissen, jannis@teunissen.net''',
        epilog='''Usage example: ./swarm_cli.py input/cs_example.txt
        -gc N2 1.0 -E_range 1e7 2e7 -E_num 10''')
    parser.add_argument('cs', type=str,
                        help='File with cross sections')
    parser.add_argument('-o', type=str, default='swarm_data.csv',
                        help='Transport data output file')
    parser.add_argument('-sigma', action='store_true',
                        help='Include standard deviations in output')
    parser.add_argument('-gc', dest='gas_comps', type=str, nargs='+',
                        required=True, metavar='gas frac',
                        help='List of gas names and fractions, '
                        'for example: N2 0.8 O2 0.2')
    parser.add_argument('-mover', type=str,
                        choices=['analytic', 'boris', 'verlet'],
                        default='analytic', help='Choice of particle mover')
    parser.add_argument('-eV_max', type=float, default=500.0,
                        help='Maximum particle energy in eV')
    parser.add_argument('-max_cpu_time', type=float, default=600.0,
                        help='Maximum CPU time per swarm')

    g = parser.add_mutually_exclusive_group()
    g.add_argument('-E', type=float,
                        help='Electric field value (V/m)')
    g.add_argument('-E_range', nargs=2, type=float, default=[1e7, 1e7],
                        help='min/max electric field (V/m)')
    parser.add_argument('-E_vary', type=str, choices=['lin', 'log'],
                        default='lin', help='How to vary electric field')
    parser.add_argument('-E_num', type=int, default=1,
                        help='Number of electric fields')

    g = parser.add_mutually_exclusive_group()
    g.add_argument('-angle', type=float,
                        help='angle between E and B (degrees)')
    g.add_argument('-angle_range', nargs=2, type=float, default=[0., 0.],
                        help='min/max angle (degrees)')
    parser.add_argument('-angle_num', type=int, default=1,
                        help='Number of angles')

    g = parser.add_mutually_exclusive_group()
    g.add_argument('-B', type=float,
                        help='Magnetic field (T)')
    g.add_argument('-B_range', nargs=2, type=float, default=[0., 0.],
                        help='min/max magnetic field (T)')
    parser.add_argument('-B_vary', type=str, choices=['lin', 'log'],
                        default='lin', help='How to vary magnetic field')
    parser.add_argument('-B_num', type=int, default=1,
                        help='Number of magnetic fields')

    parser.add_argument('-T', type=float, default=300.,
                        metavar='temperature', help='Gas temperature (K)')
    parser.add_argument('-p', type=float, default=1.0,
                        metavar='pressure', help='Gas pressure (bar)')

    parser.add_argument('-acc_v2', type=float, nargs=2, default=(1e-3, 0.0),
                        metavar=('rel', 'abs'),
                        help='Required rel./abs. error in velocity squared')
    parser.add_argument('-acc_v', type=float, nargs=2, default=(5e-3, 0.0),
                        metavar=('rel', 'abs'),
                        help='Required rel./abs. error in velocity')
    parser.add_argument('-acc_D', type=float, nargs=2, default=(1e-2, 0.0),
                        metavar=('rel', 'abs'),
                        help='Required rel./abs. error in diff. coeff')
    parser.add_argument('-acc_a', type=float, nargs=2, default=(5e-3, 10.0),
                        metavar=('rel', 'abs'),
                        help='Required rel./abs. error in alpha')
    return parser.parse_args()


def create_swarm_cfg(tmpdir, args):
    fname = tmpdir + '/base_cfg.txt'
    f = open(fname, 'w')
    f.write('output_dir = ' + tmpdir + '\n')
    f.write('gas_file = ' + args.cs + '\n')
    f.write('swarm_name = swarm_cli\n')
    f.write('gas_components = ' + ' '.join(args.gas_comps[0::2]) + '\n')
    f.write('gas_fractions = ' + ' '.join(args.gas_comps[1::2]) + '\n')
    f.write('gas_pressure = ' + str(args.p) + '\n')
    f.write('gas_temperature = ' + str(args.T) + '\n')
    f.write('acc_velocity_sq = ' + ' '.join(map(str, args.acc_v2)) + '\n')
    f.write('acc_velocity = ' + ' '.join(map(str, args.acc_v)) + '\n')
    f.write('acc_diffusion = ' + ' '.join(map(str, args.acc_D)) + '\n')
    f.write('acc_alpha = ' + ' '.join(map(str, args.acc_a)) + '\n')
    f.write('electric_field = ' + str(args.E_range[0]) + '\n')
    f.write('magnetic_field = ' + str(args.B_range[0]) + '\n')
    f.write('field_angle_degrees = ' + str(args.angle_range[0]) + '\n')
    f.write('particle_max_energy_ev = ' + str(args.eV_max) + '\n')
    f.write('max_cpu_time = ' + str(args.max_cpu_time) + '\n')
    f.write('particle_mover = ' + args.mover + '\n')
    f.close()
    return fname


def create_var_cfg(tmpdir, index, varnames, values):
    fname = tmpdir + '/var_' + str(index) + '.txt'
    f = open(fname, 'w')
    for name, val in zip(varnames, values):
        f.write(name + ' = ' + str(val) + '\n')
    f.close()
    return fname


def progress_bar(pct):
    n = int(pct * 0.4)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:40s}] {1:.1f}%'.format('=' * n, pct))
    sys.stdout.flush()


def print_swarm_info(args, E_list, B_list, angle_list):
    print("Starting particle swarm simulation")
    print("----------------------------------------")
    print("Cross section file : {}".format(args.cs))
    print("Gas components     : {}".format(' '.join(args.gas_comps[0::2])))
    print("Gas fractions      : {}".format(' '.join(args.gas_comps[1::2])))
    print("Temperature (K)    : {}".format(args.T))
    print("Pressure (bar)     : {}".format(args.p))
    print("Particle mover     : {}".format(args.mover))
    print("Output file        : {}".format(args.o))
    print("----------------------------------------")
    Evals = ["{:.2E}".format(val) for val in E_list]
    Bvals = ["{:.2E}".format(val) for val in B_list]
    Avals = ["{:.2E}".format(val) for val in angle_list]
    print("Electric fields (V/m) : {}".format(' '.join(Evals)))
    print("Magnetic fields (T)   : {}".format(' '.join(Bvals)))
    print("Angles (degrees)      : {}".format(' '.join(Avals)))
    print("----------------------------------------")


if __name__ == '__main__':
    args = get_args()

    # Set ranges if a scalar argument was given
    if args.E is not None:
        args.E_range = [args.E, args.E]
        args.E_num = 1
    if args.B is not None:
        args.B_range = [args.B, args.B]
        args.B_num = 1
    if args.angle is not None:
        args.angle_range = [args.angle, args.angle]
        args.angle_num = 1

    # Generate lists of E, angle and B values
    n_runs = args.E_num * args.angle_num * args.B_num
    if args.E_vary == 'lin':
        E_list = np.linspace(args.E_range[0], args.E_range[1], args.E_num)
    elif args.E_vary == 'log':
        E_list = np.logspace(math.log10(args.E_range[0]),
                             math.log10(args.E_range[1]), args.E_num)

    angle_list = np.linspace(args.angle_range[0], args.angle_range[1],
                             args.angle_num)

    if args.B_vary == 'lin':
        B_list = np.linspace(args.B_range[0], args.B_range[1], args.B_num)
    elif args.B_vary == 'log':
        B_list = np.logspace(math.log10(args.B_range[0]),
                             math.log10(args.B_range[1]), args.B_num)

    print_swarm_info(args, E_list, B_list, angle_list)

    try:
        tmpdir = tempfile.mkdtemp(dir=os.getcwd())
        swarm_cli_path = os.path.abspath(os.path.realpath(__file__))
        particle_swarm_exec_path = os.path.join(os.path.dirname(swarm_cli_path), "particle_swarm")

        base_cfg = create_swarm_cfg(tmpdir, args)
        cmd_list = []
        swarm_data = []

        i = 0
        names = ['electric_field', 'field_angle_degrees', 'magnetic_field']
        n_runs = len(B_list) * len(E_list) * len(angle_list)

        for B in B_list:
            for E in E_list:
                for angle in angle_list:
                    progress_bar(100. * i / n_runs)
                    i += 1
                    var_cfg = create_var_cfg(tmpdir, i, names, [E, angle, B])
                    res = check_output([particle_swarm_exec_path, base_cfg, var_cfg])
                    swarm_data.append(res)
    finally:
        progress_bar(100.)
        print("")
        shutil.rmtree(tmpdir)

    # Find transport data names
    td_names = []
    for line in swarm_data[0].splitlines():
        name = line.decode('ascii')[:40].strip()
        td_names.append(name)

    n_cols = len(td_names)
    data = np.zeros((n_runs, n_cols))
    sigma = np.zeros((n_runs, n_cols))

    for i, res in enumerate(swarm_data):
        for j, line in enumerate(res.splitlines()):
            values = line[40:].split()
            data[i, j] = values[0]
            sigma[i, j] = values[1]

    if args.sigma:
        df = pd.DataFrame(data=np.hstack([data, sigma]),
                          columns=td_names + ['stddev ' + x for x in td_names])
    else:
        df = pd.DataFrame(data=data, columns=td_names)

    df.to_csv(args.o, index=False)
    print("Results written to: ", args.o)
