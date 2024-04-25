#!/usr/bin/env python

# Author: Jannis Teunissen

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from difflib import get_close_matches

# Get and parse the command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='''Plot swarm data''')
parser.add_argument('in_file', type=argparse.FileType('r'),
                    help='CSV file with output from swarm_cli.py')
parser.add_argument('parameter', type=str, default='Mean energy (eV)',
                    help='Which quantity to plot')
parser.add_argument('-x', type=str, default=['Electric field / N (Td)'],
                    help='Use this as x-axis')
args = parser.parse_args()

df = pd.read_csv(args.in_file)

parameters = df.columns.values
x = args.x
y = args.parameter

if x not in parameters:
    x = get_close_matches(x, parameters, cutoff=0.1)[0]
    print(f'Assuming you mean "{args.x}" = "{x}"')
if y not in parameters:
    y = get_close_matches(y, parameters, cutoff=0.1)[0]
    print(f'Assuming you mean "{args.parameter}" = "{y}"')

df.set_index(x, inplace=True)

df[y].plot()
plt.legend()
plt.show()
