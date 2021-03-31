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
    g.add_argument('-EN', type=float,
                        help='Reduced electric field value (Td)')
    g.add_argument('-EN_range', nargs=2, type=float, default=[370, 370],
                        help='min/max reduced electric field (Td)')

    g = parser.add_mutually_exclusive_group()
    g.add_argument('-E_vary', type=str, choices=['lin', 'log'],
                        default='lin', help='How to vary electric field')
    g.add_argument('-EN_vary', type=str, choices=['lin', 'log'],
                        default='lin', help='How to vary reduced electric field')

    g = parser.add_mutually_exclusive_group()
    g.add_argument('-E_num', type=int, default=1,
                        help='Number of electric fields')
    g.add_argument('-EN_num', type=int, default=1,
                        help='Number of reduced electric fields')

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

    parser.add_argument('-base_cfg', type=str, help='Optional configuration file which fixes variables across runs (e.g. when using E_range).')
    
    return parser.parse_args()


def read_cfg_to_dict(cfg_filename):
    cfg_dict = {}
    with open(cfg_filename) as f:
        for line in f.readlines():
            sanitized_string = line.lstrip()
            if len(sanitized_string) > 0:
                if sanitized_string[0] != "#":
                    split_line = sanitized_string.split("=")
                    key = split_line[0].strip()

                    val = split_line[1].strip()
                    if len(val.split()) > 1:
                        val = [x for x in val.split()]

                    cfg_dict[key] = val

    return cfg_dict

def create_swarm_cfg(tmpdir, args):
    fname = tmpdir + '/base_cfg.txt'

    if len(args.base_cfg) > 0:
        cfg_dict = read_cfg_to_dict(args.base_cfg)

        # Write to a temporary cfg file any config variable NOT in the passed base cfg
        # Note: Since gas components and gas fractions are mandatory in the cli
        # the ones written in the passed base cfg will NOT be used. The ones
        # from the cli will be used and will be written to a temporary cfg file
        with open(fname, "w") as f:
            f.write(f"output_dir = {tmpdir}\n")
            f.write(f"gas_file = {args.cs}\n")

            if "swarm_name" not in cfg_dict:
                f.write("swarm_name = swarm_cli\n")
            
            f.write(f"gas_components = {' '.join(args.gas_comps[0::2])}\n")
            f.write(f"gas_fractions = {' '.join(args.gas_comps[1::2])}\n")

            if "gas_pressure" not in cfg_dict:
                f.write(f"gas_pressure = {args.p}\n")
            if "gas_temperature" not in cfg_dict:
                f.write(f"gas_temperature = {args.T}\n")
            if "acc_velocity_sq" not in cfg_dict:
                f.write(f"acc_velocity_sq = {' '.join(map(str, args.acc_v2))}\n")
            if "acc_velocity" not in cfg_dict:
                f.write(f"acc_velocity = {' '.join(map(str, args.acc_v))}\n")
            if "acc_diffusion" not in cfg_dict:
                f.write(f"acc_diffusion = {' '.join(map(str, args.acc_D))}\n")
            if "acc_alpha" not in cfg_dict:
                f.write(f"acc_alpha = {' '.join(map(str, args.acc_a))}\n")
            if "electric_field" not in cfg_dict:
                f.write(f"electric_field = {args.E_range[0]}\n")
            if "magnetic_field" not in cfg_dict:
                f.write(f"magnetic_field = {args.B_range[0]}\n")
            if "field_angle_degrees" not in cfg_dict:
                f.write(f"field_angle_degrees = {args.angle_range[0]}\n")
            if "particle_max_energy_ev" not in cfg_dict:
                f.write(f"particle_max_energy_ev = {args.eV_max}\n")
            if "max_cpu_time" not in cfg_dict:
                f.write(f"max_cpu_time = {args.max_cpu_time}\n")
            if "particle_mover" not in cfg_dict:
                f.write(f"particle_mover = {args.mover}\n")

        return [args.base_cfg, fname]
    else:
        with open(fname, "w") as f:
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
    if args.base_cfg is not None:
        cfg_dict = read_cfg_to_dict(args.base_cfg)
        if "gas_temperature" in cfg_dict:
            args.T = cfg_dict["gas_temperature"]
        if "gas_pressure" in cfg_dict:
            args.p = cfg_dict["gas_pressure"]
        if "particle_mover" in cfg_dict:
            args.mover = cfg_dict["particle_mover"]

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

# Convert EN in Td to E in V / m
def EN_to_E(EN, pressure_bar):
    Td = 1e-21  # V m2
    n0 = 2.414321e25  # m-3

    return EN * Td * (n0 * pressure_bar)

if __name__ == '__main__':
    args = get_args()

    pressure_bar = args.p
    if args.base_cfg is not None:
        cfg_dict = read_cfg_to_dict(args.base_cfg)
        pressure_bar = float(cfg_dict["gas_pressure"])

    # Set E_list
    if args.E is not None:
        E_list = np.array([args.E])
    elif args.E_num > 1:
        if args.E_vary == 'lin':
            E_list = np.linspace(args.E_range[0], args.E_range[1], args.E_num)
        elif args.E_vary == 'log':
            E_list = np.logspace(math.log10(args.E_range[0]),
                                math.log10(args.E_range[1]), args.E_num)
    elif args.EN is not None:
        E_list = np.array([EN_to_E(args.EN, pressure_bar)])
    elif args.EN_num > 1:
        if args.EN_vary == 'lin':
            E_list = np.linspace(EN_to_E(args.EN_range[0], pressure_bar), EN_to_E(args.EN_range[1], pressure_bar), args.EN_num)
        elif args.EN_vary == 'log':
            E_list = np.logspace(math.log10(EN_to_E(args.EN_range[0], pressure_bar)),
                                math.log10(EN_to_E(args.EN_range[1], pressure_bar)), args.EN_num)
    else:
        # When nothing has been passed for electric field
        E_list = np.array([args.E_range[0]])
    
    # Set B_list
    if args.B is not None:
        B_list = np.array([args.B])
    elif args.B_num > 1:
        if args.B_vary == 'lin':
            B_list = np.linspace(args.B_range[0], args.B_range[1], args.B_num)
        elif args.B_vary == 'log':
            B_list = np.logspace(math.log10(args.B_range[0]),
                                math.log10(args.B_range[1]), args.B_num)
    else:
        # When nothing has been passed for B field
        B_list = np.array([args.B_range[0]])

    # Set angle_list
    if args.angle is not None:
        angle_list = np.array([args.angle])
    elif args.angle_num > 1:
        angle_list = np.linspace(args.angle_range[0], args.angle_range[1],
                                 args.angle_num)
    else:
        # When nothing has been passed for angles
        angle_list = np.array([args.angle_range[0]])


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
                    if type(base_cfg) == list:
                        res = check_output([particle_swarm_exec_path, *base_cfg, var_cfg])
                    else:
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
