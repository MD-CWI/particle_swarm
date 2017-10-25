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
import time
import sys
import math
import numpy as np
from multiprocessing import Pool, Manager, cpu_count
from subprocess import check_output, CalledProcessError


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
    parser.add_argument('-of', type=str, default='results.txt',
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
                        help='min/max angle (degrees)')
    g.add_argument('-angle_range', nargs=2, type=float, default=[0., 0.],
                        help='min/max angle (degrees)')
    parser.add_argument('-angle_vary', type=str, choices=['lin', 'log'],
                        default='lin', help='How to vary angle')
    parser.add_argument('-angle_num', type=int, default=1,
                        help='Number of angles')

    g = parser.add_mutually_exclusive_group()
    g.add_argument('-B', type=float,
                        help='Magnetic field range (T)')
    g.add_argument('-B_range', nargs=2, type=float, default=[0., 0.],
                        help='Magnetic field range (T)')
    parser.add_argument('-B_vary', type=str, choices=['lin', 'log'],
                        default='lin', help='How to vary magnetic field')
    parser.add_argument('-B_num', type=int, default=1,
                        help='Number of magnetic fields')

    parser.add_argument('-T', type=float, default=300.,
                        metavar='temperature', help='Gas temperature (K)')
    parser.add_argument('-p', type=float, default=1.0,
                        metavar='pressure', help='Gas pressure (bar)')
    parser.add_argument('-np', type=int, default=max(1, cpu_count()//2),
                        help='Number of parallel proccesses')

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
    f.write('consecutive_run = T\n')
    f.write('dry_run = F\n')
    f.write('acc_velocity_sq = ' + ' '.join(map(str, args.acc_v2)) + '\n')
    f.write('acc_velocity = ' + ' '.join(map(str, args.acc_v)) + '\n')
    f.write('acc_diffusion = ' + ' '.join(map(str, args.acc_D)) + '\n')
    f.write('acc_alpha = ' + ' '.join(map(str, args.acc_a)) + '\n')
    f.write('electric_field = ' + str(args.E_range[0]) + '\n')
    f.write('magnetic_field = ' + str(args.B_range[0]) + '\n')
    f.write('field_angle_degrees = ' + str(args.angle_range[0]) + '\n')
    f.write('particle_mover = ' + args.mover + '\n')
    f.close()
    return fname


def create_init_cfg(tmpdir):
    fname = tmpdir + '/init_cfg.txt'
    f = open(fname, 'w')
    f.write('consecutive_run = F\n')
    f.write('dry_run = T\n')
    f.close()
    return fname


def create_var_cfg(tmpdir, index, varnames, values):
    fname = tmpdir + '/var_' + str(index) + '.txt'
    f = open(fname, 'w')
    for name, val in zip(varnames, values):
        f.write(name + ' = ' + str(val) + '\n')
    f.close()
    return fname


def pswarm_wrapper(cmd_and_num):
    cmd, num = cmd_and_num
    try:
        res = check_output(cmd)
    except:
        print("Error when running: " + cmd)
        sys.exit(1)
    num.value += 1
    return res


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
    print("Output file        : {}".format(args.of))
    print("Number of CPUs     : {}".format(args.np))
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

    if args.angle_vary == 'lin':
        angle_list = np.linspace(args.angle_range[0], args.angle_range[1],
                                 args.angle_num)
    elif args.angle_vary == 'log':
        angle_list = np.logspace(math.log10(args.angle_range[0]),
                             math.log10(args.angle_range[1]), args.angle_num)

    if args.B_vary == 'lin':
        B_list = np.linspace(args.B_range[0], args.B_range[1], args.B_num)
    elif args.B_vary == 'log':
        B_list = np.logspace(math.log10(args.B_range[0]),
                             math.log10(args.B_range[1]), args.B_num)

    print_swarm_info(args, E_list, B_list, angle_list)

    try:
        tmpdir = tempfile.mkdtemp(dir=os.getcwd())

        base_cfg = create_swarm_cfg(tmpdir, args)
        init_cfg = create_init_cfg(tmpdir)

        # Perform a dry run of particle_swarm to generate lookup tables for the
        # cross sections
        try:
            out = check_output(['./particle_swarm', base_cfg, init_cfg])
        except CalledProcessError, e:
            print("particle_swarm returned an error")
            print(e.output)
            sys.exit(1)

        with open(tmpdir + '/swarm_cli_cs_summary.txt') as f:
            cs_info = f.read()
            print(cs_info.rstrip())
            print("------------------------------------------------------")

        cmd_list = []

        pool = Pool(processes=args.np)
        mgr = Manager()
        num = mgr.Value('i', 0)
        i = 0
        names = ['electric_field', 'field_angle_degrees', 'magnetic_field']

        for E in E_list:
            for angle in angle_list:
                for B in B_list:
                    i += 1
                    var_cfg = create_var_cfg(tmpdir, i, names, [E, angle, B])
                    cmd_list.append([['./particle_swarm',
                                      base_cfg, var_cfg], num])

        res = pool.map_async(pswarm_wrapper, cmd_list)
        while num.value < n_runs:
            time.sleep(1.)
            progress_bar(num.value / n_runs * 100)
        print('')
        swarm_data = res.get()

    finally:
        shutil.rmtree(tmpdir)

    # Find transport data names
    td_names = []
    for line in swarm_data[0].splitlines():
        name = line.split()[0].decode('ascii')
        td_names.append(name)

    n_cols = len(td_names)

    if args.sigma:
        td_matrix = np.zeros((n_runs, 2 * n_cols))
    else:
        td_matrix = np.zeros((n_runs, n_cols))

    for i, res in enumerate(swarm_data):
        for j, line in enumerate(res.splitlines()):
            # Format is name = value stddev
            # line.split()[2] corresponds to 'value'
            # line.split()[3] corresponds to 'stddev'
            td_matrix[i, j] = line.split()[2]
            if args.sigma:
                td_matrix[i, j + n_cols] = line.split()[3]

    header = '# ' + ' '.join(td_names)

    # Prepend 's' to indicated standard deviations
    if args.sigma:
        snames = ['s' + name for name in td_names]
        header += ' ' + ' '.join(snames)

    # Write header manually to support numpy < 1.7
    with open(args.of, 'wb') as f:
        f.write(header.encode('ascii'))
        f.write('\n'.encode('ascii'))
        np.savetxt(f, td_matrix, fmt=b'%10.3e')

    print("Results written to: ", args.of)
