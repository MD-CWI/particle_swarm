#!/usr/bin/env python

# Author: Jannis Teunissen

import argparse
import pandas as pd
import matplotlib.pyplot as plt

parameters = [
    'Gas number density (1/m3)',
    'Electric field (V/m)',
    'Electric field / N (Td)',
    'Mean energy (eV)',
    'Townsend ioniz. coef. alpha/N (m2)',
    'Townsend attach. coef. eta/N (m2)',
    'Mobility (m2/V/s)',
    'Mobility *N (1/m/V/s)',
    'Bulk mobility (m2/V/s)',
    'Bulk mobility *N (1/m/V/s)',
    'Flux L diffusion coef. (m2/s)',
    'Flux L diffusion coef. *N (1/m/s)',
    'Bulk L diffusion coef. (m2/s)',
    'Bulk L diffusion coef. *N (1/m/s)',
    'Flux T diffusion coef. (m2/s)',
    'Flux T diffusion coef. *N (1/m/s)',
    'Bulk T diffusion coef. (m2/s)',
    'Bulk T diffusion coef. *N (1/m/s)',
    'Townsend ioniz. coef. alpha (1/m)',
    'Townsend attach. coef. eta (1/m)',
    'Total ionization freq. (1/s)',
    'Total attachment freq. (1/s)',
    'Total collision rate (1/s)',
    'Drift velocity x (m/s)',
    'Drift velocity y (m/s)',
    'Drift velocity z (m/s)',
    'Bulk drift velocity x (m/s)',
    'Bulk drift velocity y (m/s)',
    'Bulk drift velocity z (m/s)',
    'Diffusion coefficient x (m2/s)',
    'Diffusion coefficient y (m2/s)',
    'Diffusion coefficient z (m2/s)',
    'Bulk diffusion coef. x (m2/s)',
    'Bulk diffusion coef. y (m2/s)',
    'Bulk diffusion coef. z (m2/s)',
    'Velocity squared x (m2/s2)',
    'Velocity squared y (m2/s2)',
    'Velocity squared z (m2/s2)'
]

# Get and parse the command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='''Plot swarm data''')
parser.add_argument('in_file', type=argparse.FileType('r'),
                    help='CSV file with output from swarm_cli.py')
parser.add_argument('parameter', type=str, default='Mean energy (eV)',
                    choices=parameters, help='Which quantity to plot')
parser.add_argument('-x', type=str, default=['Electric field / N (Td)'],
                    choices=parameters, help='Use this as x-axis')
args = parser.parse_args()

df = pd.read_csv(args.in_file)
df.set_index(args.x, inplace=True)

df[args.parameter].plot()
plt.legend()
plt.show()
