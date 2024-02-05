#!/usr/bin/env python3

# Plot electron energy relaxation time computed from mobility
# Author: Jannis Teunissen

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Get and parse the command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='''Plot relaxation times''')
parser.add_argument('in_file', type=argparse.FileType('r'),
                    help='CSV file with output from swarm_cli.py')
parser.add_argument('-savefig', type=str,
                    help='Save figure into file')
parser.add_argument('-ymax', type=float,
                    help='Max. of y scale')
args = parser.parse_args()

df = pd.read_csv(args.in_file)

fig, ax = plt.subplots(figsize=(6, 4), layout='constrained')

x = 'E (MV/m)'
df[x] = df['Electric field (V/m)'] * 1e-6
df['energy gain (eV/s)'] = df['Mobility (m2/V/s)'] * \
    df['Electric field (V/m)']**2
df['tau energy'] = df['Mean energy (eV)'] / df['energy gain (eV/s)']
df['tau coll'] = 1 / df['Total collision rate (1/s)']
df['tau ionization'] = 1 / (np.abs(df['Townsend ioniz. coef. alpha (1/m)'] -
                                   df['Townsend attach. coef. eta (1/m)']) *
                            df['Mobility (m2/V/s)'] *
                            df['Electric field (V/m)'])

ax.plot(df[x], df['tau energy'], label=r'$\tau_\varepsilon$')
ax.plot(df[x], df['tau coll'], label=r'$\tau_\mathrm{coll}$')
ax.plot(df[x], df['tau ionization'], label=r'$\tau_\alpha$')

ax.set_xlabel(x)
ax.set_ylabel('t (s)')
ax.legend()
ax.set_yscale('log')

if args.ymax:
    ax.set_ylim(None, args.ymax)

if args.savefig:
    plt.savefig(args.savefig, dpi=300, bbox_inches='tight')
else:
    plt.show()
