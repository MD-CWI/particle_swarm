#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import numpy as np


def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('in_file', type=argparse.FileType('r'),
                        help='File with output from swarm_cli')
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    # Read header
    hdr = args.in_file.readline()

    # Strip comment and newline
    hdr = hdr[2:].strip()
    colnames = hdr.split(' ')
    print(colnames)

    # Read data
    td = np.loadtxt(args.in_file)

