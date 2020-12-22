#!/usr/bin/env python

# Author: Jannis Teunissen

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import pandas as pd
import datetime
import re


def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Convert particle swarm data to tables that can be used in discharge
        simulations with the local field approximation. Author: Jannis
        Teunissen, jannis@teunissen.net''',
        epilog='Example: ./swarm_to_lfa.py swarm_data.csv AIR > table.txt')
    parser.add_argument('in_file', type=argparse.FileType('r'),
                        help='CSV file with output from swarm_cli.py')
    parser.add_argument('-output_format', type=str, default='afivo-streamer',
                        choices=['streamer_1d', 'afivo-streamer'],
                        help='Generate transport data for this program')
    return parser.parse_args()


def write_col(df, x, y):
    print(y)
    print("-----------------------")
    for x, y in zip(df[x], df[y]):
        print("{:.8E} {:.8E}".format(x, y))
    print("-----------------------")
    print("")


if __name__ == '__main__':
    args = get_args()

    # Read data
    df = pd.read_csv(args.in_file)

    print('# Generated by swarm_to_lfa.py')
    print('# Date:', datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))

    if args.output_format == 'afivo-streamer':
        x = 'Electric field / N (Td)'
        print('# First column: ', x, '\n')

        print('# Guess for reaction list')
        print('reaction_list')
        print('-----------------------')
        print('# Ionization')
        pattern = re.compile(r'^C\d+ (\S+) Ionization.*$')
        for c in df.columns.values:
            m = pattern.match(c)
            if m:
                print('e + ' + m.group(1) + ' -> e + e + ' +
                      m.group(1) + '+,field_table,' + m.group(0))
        print('# Attachment')
        pattern = re.compile(r'^C\d+ (\S+) Attachment.*$')
        for c in df.columns.values:
            m = pattern.match(c)
            if m:
                print('e + ' + m.group(1) + ' -> ' +
                      m.group(1) + '-,field_table,' + m.group(0))
        print('-----------------------\n')

        write_col(df, x, 'Mean energy (eV)')
        write_col(df, x, 'Mobility *N (1/m/V/s)')
        write_col(df, x, 'Bulk mobility *N (1/m/V/s)')
        write_col(df, x, 'Townsend ioniz. coef. alpha/N (m2)')
        write_col(df, x, 'Townsend attach. coef. eta/N (m2)')
        write_col(df, x, 'Flux L diffusion coef. *N (1/m/s)')
        write_col(df, x, 'Bulk L diffusion coef. *N (1/m/s)')
        write_col(df, x, 'Flux T diffusion coef. *N (1/m/s)')
        write_col(df, x, 'Bulk T diffusion coef. *N (1/m/s)')
        write_col(df, x, 'Total ionization freq. (1/s)')
        write_col(df, x, 'Total attachment freq. (1/s)')

        # Print rates for other collisions
        pattern = re.compile(r'^C\d+')
        for c in df.columns.values:
            if pattern.match(c):
                write_col(df, x, c)

    elif args.output_format == 'streamer_1d':
        x = 'Electric field (V/m)'
        print('# First column: ', x, '\n')

        write_col(df, x, 'Mobility (m2/V/s)')
        write_col(df, x, 'Bulk mobility (m2/V/s)')
        write_col(df, x, 'Townsend ioniz. coef. alpha (1/m)')
        write_col(df, x, 'Townsend attach. coef. eta (1/m)')
        write_col(df, x, 'Mean energy (eV)')
        write_col(df, x, 'Flux L diffusion coef. (m2/s)')
        write_col(df, x, 'Bulk L diffusion coef. (m2/s)')
        write_col(df, x, 'Flux T diffusion coef. (m2/s)')
        write_col(df, x, 'Bulk T diffusion coef. (m2/s)')
