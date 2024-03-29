#!/usr/bin/env python

# Script to plot cross section data

import numpy as np
import argparse
import scipy.constants as const
from bolos_parser import parse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Plot cross sections')
parser.add_argument('cs', type=str, help='File with cross sections')
parser.add_argument('-collision_rate', action='store_true',
                    help='Show collision rates instead')
parser.add_argument('-only_show', type=int, nargs='+',
                    help='Index of cross sections to show')
parser.add_argument('-xlim', type=float, nargs=2,
                    help='Min/max energy')
parser.add_argument('-yscale', type=str, default='linear',
                    help='Type of y-scale')
parser.add_argument('-savefig', type=str,
                    help='File name for figure')
args = parser.parse_args()

with open(args.cs, 'r') as f:
    all_cs = parse(f)

if args.only_show:
    all_cs = [all_cs[i] for i in args.only_show]

fig, ax = plt.subplots(figsize=(5, 5), layout='constrained')

for cs in all_cs:
    cs_data = np.array(cs['data'])
    x = cs_data[:, 0]
    y = cs_data[:, 1]

    label = cs['kind'].lower()
    if 'product' in cs:
        label += ' ' + cs['product']

    if args.collision_rate:
        y = np.sqrt(2 * x * const.electron_volt/const.electron_mass) * y
    ax.plot(x, y, label=label)

ax.set_xlabel('Energy (eV)')

if args.xlim:
    ax.set_xlim(args.xlim)

if args.collision_rate:
    ax.set_ylabel('Collision rate (s$^{-1}$ N$^{-1}$)')
else:
    ax.set_ylabel('Cross section (m$^2$)')

plt.legend()

ax.set_yscale(args.yscale)

if args.savefig:
    plt.savefig(args.savefig, bbox_inches='tight', dpi=300)
else:
    plt.show()
