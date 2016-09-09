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
from subprocess import check_output


def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        description='''Command line interface to compute electron
        transport data from swarm simulations.
        Author: Jannis Teunissen, jannis@teunissen.net''',
        epilog='''Example:
        ./swarm_cli.py crosssec.txt -gc N2 1.0 -vary E -vlin 1e7 2e7 10''')
    parser.add_argument('cs', type=str,
                        help='File with cross sections')
    parser.add_argument('-out_file', type=str, default='results.txt',
                        help='Output file, default is results.txt')
    parser.add_argument('-gc', dest='gas_comps', type=str, nargs='+',
                        required=True, metavar='gas frac',
                        help='List of gas names and fractions, '
                        'for example: N2 0.8 O2 0.2')
    parser.add_argument('-vary', type=str, choices=['E', 'B', 'angle'],
                        default='E', required=True,
                        help='Which quantity to vary')
    parser.add_argument('-mover', type=str,
                        choices=['analytic', 'boris', 'verlet'],
                        default='verlet', help='Choice of particle mover')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-vlist', type=float, nargs='+', metavar='val',
                       help='List of values to use')
    group.add_argument('-vlin', type=float, nargs=3,
                       metavar=('min', 'max', 'N'),
                       help='Linear range of values to use')
    group.add_argument('-vlog', type=float, nargs=3,
                       metavar=('min', 'max', 'N'),
                       help='Logarithmic range of values to use')

    parser.add_argument('-E', type=float, default=1e7,
                         help='Electric field (V/m)')
    parser.add_argument('-B', type=float, default=0.0,
                         help='Magnetic field (V/m)')
    parser.add_argument('-angle', type=float, default=0.0,
                         help='Angle between E and B (degrees)')

    parser.add_argument('-T', type=float, default=300.,
                        metavar='temperature', help='Gas temperature (K)')
    parser.add_argument('-p', type=float, default=1.0,
                        metavar='pressure', help='Gas pressure (bar)')
    parser.add_argument('-np', type=int, default=cpu_count(),
                        help='Number of parallel proccesses')

    parser.add_argument('-en', type=float, nargs=2, default=(1e-3, 0.0),
                        metavar=('rel', 'abs'),
                        help='Required rel./abs. error in energy')
    parser.add_argument('-mu', type=float, nargs=2, default=(5e-3, 0.0),
                        metavar=('rel', 'abs'),
                        help='Required rel./abs. error in mobility')
    parser.add_argument('-D', type=float, nargs=2, default=(1e-2, 0.0),
                        metavar=('rel', 'abs'),
                        help='Required rel./abs. error in diff. coeff')
    parser.add_argument('-a', type=float, nargs=2, default=(5e-3, 10.0),
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
    f.write('acc_energy = ' + ' '.join(map(str, args.en)) + '\n')
    f.write('acc_mobility = ' + ' '.join(map(str, args.mu)) + '\n')
    f.write('acc_diffusion = ' + ' '.join(map(str, args.D)) + '\n')
    f.write('acc_alpha = ' + ' '.join(map(str, args.a)) + '\n')
    f.write('electric_field = ' + str(args.E) + '\n')
    f.write('magnetic_field = ' + str(args.B) + '\n')
    f.write('field_angle_degrees = ' + str(args.angle) + '\n')
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


def create_var_cfg(tmpdir, index, varname, value):
    fname = tmpdir + '/var_' + str(index) + '.txt'
    f = open(fname, 'w')
    f.write(varname + ' = ' + str(value) + '\n')
    f.close()
    return fname


def pswarm_wrapper(cmd_and_num):
    cmd, num = cmd_and_num
    res = check_output(cmd)
    num.value += 1
    return res


def progress_bar(pct):
    n = int(pct * 0.4)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:40s}] {1:.1f}%'.format('=' * n, pct))
    sys.stdout.flush()


if __name__ == '__main__':
    args = get_args()

    # Make sure we have a list of electric fields
    if args.vlin:
        args.vlist = np.linspace(args.vlin[0], args.vlin[1], args.vlin[2])
    elif args.vlog:
        args.vlist = np.logspace(math.log10(args.vlog[0]),
                                 math.log10(args.vlog[1]), args.vlog[2])
    n_flds = len(args.vlist)

    try:
        tmpdir = tempfile.mkdtemp(dir=os.getcwd())

        base_cfg = create_swarm_cfg(tmpdir, args)
        init_cfg = create_init_cfg(tmpdir)

        # Perform a dry run of particle_swarm to generate lookup tables for the
        # cross sections
        check_output(['./particle_swarm', base_cfg, init_cfg])
        cmd_list = []

        pool = Pool(processes=args.np)
        mgr = Manager()
        num = mgr.Value('i', 0)

        for i, value in enumerate(args.vlist):
            if args.vary == 'E':
                var_cfg = create_var_cfg(tmpdir, i, 'electric_field', value)
            elif args.vary == 'B':
                var_cfg = create_var_cfg(tmpdir, i, 'magnetic_field', value)
            elif args.vary == 'angle':
                var_cfg = create_var_cfg(tmpdir, i, 'field_angle_degrees', value)
            cmd_list.append([['./particle_swarm', base_cfg, var_cfg], num])

        res = pool.map_async(pswarm_wrapper, cmd_list)
        while num.value < n_flds:
            time.sleep(1.)
            progress_bar(num.value / n_flds * 100)
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
    td_matrix = np.zeros((n_flds, n_cols))

    for i, res in enumerate(swarm_data):
        for j, line in enumerate(res.splitlines()):
            td_matrix[i, j] = line.split()[2]

    header = ' '.join(td_names)
    np.savetxt(args.out_file, td_matrix, fmt=b'%10.3e', header=header)
    print("Done, results written to: ", args.out_file)
