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
parser.add_argument('-N0', type=float, default=1,
                    help='Gas number density (1/m^3)')
parser.add_argument('-only_show', type=int, nargs='+',
                    help='Index of cross sections to show')
parser.add_argument('-show_sum', action='store_true',
                    help='Show sum of the cross sections')
parser.add_argument('-xlim', type=float, nargs=2,
                    help='Min/max energy')
parser.add_argument('-show_legend', action='store_true',
                    help='Show legend')
parser.add_argument('-yscale', type=str, default='linear',
                    help='Type of y-scale')
parser.add_argument('-xscale', type=str, default='linear',
                    help='Type of x-scale')
parser.add_argument('-savefig', type=str,
                    help='File name for figure')
args = parser.parse_args()

with open(args.cs, 'r') as f:
    all_cs = parse(f)

if args.only_show:
    all_cs = [all_cs[i] for i in args.only_show]

if args.show_sum:
    if args.xlim is None:
        raise ValueError('Show sum requires xlim argument')
    # Quadratic spacing
    x_plot = args.xlim[0] + np.linspace(0., 1., 1000)**2 * \
        (args.xlim[1] - args.xlim[0])
    y_sum = np.zeros(1000)

fig, ax = plt.subplots(figsize=(5, 3), layout='constrained')
cs_fill_color = {'ionization': 'red',
                 'excitation': 'yellow',
                 'elastic': 'green'}

for cs in all_cs:
    cs_data = np.array(cs['data'])
    x = cs_data[:, 0]
    y = cs_data[:, 1]

    if args.collision_rate:
        # Note that this linearly interpolates between collision rates, not
        # cross sections, which is usually not what should happen
        v = np.sqrt(2 * x * const.electron_volt/const.electron_mass)
        y = v * y * args.N0

    if args.show_sum:
        y_sum_prev = y_sum.copy()
        y_sum += np.interp(x_plot, x, y, left=0., right=0.)
        y = y_sum
        x = x_plot

    cs_kind = cs['kind'].lower()
    label = cs_kind

    if 'product' in cs:
        label += ' ' + cs['product']

    ax.plot(x, y, label=label)

    if args.show_sum:
        if cs_kind in cs_fill_color:
            ax.fill_between(x, y, y_sum_prev, color=cs_fill_color[cs_kind],
                            alpha=0.5)
        else:
            p = ax.fill_between(x, y, y_sum_prev, alpha=0.5)

ax.set_xlabel('Energy (eV)')

if args.xlim:
    ax.set_xlim(args.xlim)

if args.collision_rate:
    if args.N0 != 1.0:
        ax.set_ylabel('Collision rate (s$^{-1}$)')
    else:
        ax.set_ylabel('Collision rate (s$^{-1}$ N$^{-1}$)')
else:
    ax.set_ylabel('Cross section (m$^2$)')

if args.show_legend:
    plt.legend()

ax.set_yscale(args.yscale)
ax.set_xscale(args.xscale)

if args.savefig:
    plt.savefig(args.savefig, bbox_inches='tight', dpi=300)
else:
    plt.show()
