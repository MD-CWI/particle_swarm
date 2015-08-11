Electron swarm simulator
====

This code can be used to perform electron swarm simulations, using a Monte Carlo
particle model. The Monte Carlo particle model is written in Fortran 90, and
uses some Fortran 2003 features. There is a Python 2/3 command line interface
which can be used to perform parallel swarm computations, see the example below.

### Requirements

1. Particle model -- gfortran 4.8 or newer
2. Command line interface -- Python 2.7 or newer

### Getting the code

    git clone https://github.com/jannisteunissen/particle_swarm
    cd particle_swarm

### Compiling the code

    make

### Using the command line interface

    # To see command line options
    ./swarm_cli.py --help

    # Example of invocation
    ./swarm_cli.py crosssec.txt results.txt -gc N2 1.0 -flin 1e7 2e7 10

The last command uses cross sections from the file crosssec.txt, stores results
in results.txt and performs swarm simulations in nitrogen (N2) for 10 electric
fields between 1e7 V/m and 2e7 V/m.

### Author
Jannis Teunissen, jannis@teunissen.net
