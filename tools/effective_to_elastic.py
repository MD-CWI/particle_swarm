#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from parser import parse
from scipy.interpolate import interp1d
import argparse

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('input_files', type=str, nargs='+',
               help='Input files')
p.add_argument('-savefig', type=str,
               help='Save figure with this name')
p.add_argument('-nx', type=int, default=200,
               help='Number of points for interpolation')
p.add_argument('-max_eV', type=float, default=500.,
               help='Maximal energy in eV')
p.add_argument('-spacing_order', type=float, default=2.0,
               help='Polynomial order for energy spacing (1.0: linear)')
p.add_argument('-saveconverted', type=str,
               help='Save converted effective cross section in this file')

args = p.parse_args()

cs_data = []

for fname in args.input_files:
    with open(fname, 'r') as f:
        data = parse(f)
        cs_data.append(data)

x = np.linspace(0., 1., args.nx)**(args.spacing_order) * args.max_eV

fig, axs = plt.subplots(3, figsize=(8, 12), sharex=True)

for fname, data in zip(args.input_files, cs_data):
    sum_inelastic = np.zeros(args.nx)
    total_cs = np.zeros(args.nx)
    elastic = None
    effective = None

    for d in data:
        xy = np.array(d['data'])
        f = interp1d(xy[:, 0], xy[:, 1], kind='linear',
                     bounds_error=False, fill_value=0.)
        total_cs = total_cs + f(x)

        if d['kind'] == 'ELASTIC':
            if elastic is not None:
                raise ValueError('Multiple ELASTIC in '.format(fname))
            else:
                elastic = f(x)
        elif d['kind'] == 'EFFECTIVE':
            if effective is not None:
                raise ValueError('Multiple EFFECTIVE in '.format(fname))
            else:
                effective = f(x)
        else:
            sum_inelastic = sum_inelastic + f(x)

    if (elastic is not None) == (effective is not None):
        raise ValueError(
            'Need either EFFECTIVE or ELASTIC in {}'.format(fname))

    if args.saveconverted is not None and effective is not None:
        # Find effective
        old_cs = [d for d in data if d['kind'] == 'EFFECTIVE'][0]

        with open(args.saveconverted, 'w') as f:
            f.write('ELASTIC\n')
            f.write(old_cs['target'] + '\n')
            f.write(str(old_cs['mass_ratio']) + '\n')
            f.write(old_cs['comment'] + '\n')
            f.write('COMMENT: converted from EFFECTIVE\n')
            f.write('---------------\n')
            for i in range(args.nx):
                f.write('{:<20.12} {:<20.12}'.format(
                    x[i], effective[i] - sum_inelastic[i]).strip() + '\n')
            f.write('---------------\n')
        print('Converted EFFECTIVE from {} to {}'.format(
            fname, args.saveconverted))

    if args.savefig is not None:
        if effective is not None:
            axs[0].plot(x, effective - sum_inelastic, ':',
                        label='elastic* ' + fname)
            axs[2].plot(x, effective, ':', label='effective ' + fname)
        else:
            axs[0].plot(x, elastic, label=d['kind'] + ' ' + fname)
            axs[2].plot(x, total_cs, '--', label='total ' + fname)
        axs[1].plot(x, sum_inelastic, '--', label='sum inelastic ' + fname)

if args.savefig is not None:
    for i in range(len(axs)):
        axs[i].set_xscale('log')
        axs[i].set_ylabel('cross section (m^2)')
        axs[i].set_xlabel('energy (eV)')
        axs[i].legend()

    plt.tight_layout()
    plt.savefig(args.savefig, dpi=200, bbox_inches='tight')
    print('Saved ' + args.savefig)
