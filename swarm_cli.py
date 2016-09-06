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
        ./swarm_cli.py -cs crosssec.txt -of results.txt -gc N2 1.0 -flist 1e7 2e7''')
    parser.add_argument('cs', type=str,
                        help='File with cross sections')
    parser.add_argument('-out_file', type=str, default='results.txt',
                        help='Output file, default is results.txt')
    parser.add_argument('-gc', dest='gas_comps', type=str, nargs='+',
                        required=True, metavar='gas frac',
                        help='List of gas names and fractions, '
                        'for example: N2 0.8 O2 0.2')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-flist', type=float, nargs='+', metavar='E',
                       help='List of electric fields (V/m)')
    group.add_argument('-flin', type=float, nargs=3,
                       metavar=('min', 'max', 'N'),
                       help='Linear range of N electric fields (V/m)')
    group.add_argument('-flog', type=float, nargs=3,
                       metavar=('min', 'max', 'N'),
                       help='Logarithmic range of N electric fields (V/m)')

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
    f.close()
    return fname


def create_init_cfg(tmpdir):
    fname = tmpdir + '/init_cfg.txt'
    f = open(fname, 'w')
    f.write('consecutive_run = F\n')
    f.write('dry_run = T\n')
    f.close()
    return fname


def create_fld_cfg(tmpdir, index, fld):
    fname = tmpdir + '/fld_' + str(index) + '.txt'
    f = open(fname, 'w')
    f.write('electric_field = ' + str(fld) + '\n')
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
    if args.flin:
        args.flist = np.linspace(args.flin[0], args.flin[1], args.flin[2])
    elif args.flog:
        args.flist = np.logspace(math.log10(args.flog[0]),
                                 math.log10(args.flog[1]), args.flog[2])
    n_flds = len(args.flist)

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

        for i, fld in enumerate(args.flist):
            fld_cfg = create_fld_cfg(tmpdir, i, fld)
            cmd_list.append([['./particle_swarm', base_cfg, fld_cfg], num])

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
            td_matrix[i, j] = line.split()[1]

    header = ' '.join(td_names)
    np.savetxt(args.out_file, td_matrix, fmt=b'%10.3e', header=header)
    print("Done, results written to: ", args.out_file)
