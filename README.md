Electron swarm simulator
====

This code can be used to perform electron swarm simulations, using a Monte Carlo
particle model. The Monte Carlo particle model is written in Fortran 90, and
uses some Fortran 2003 features. There is a Python 2/3 command line interface
which can be used to perform parallel swarm computations, see the example below.

### Installation & compilation

Requirements
Particle model: gfortran 4.8 or newer
Command line interface: python 2.7 or newer

Getting the code and compiling

    git clone https://github.com/jannisteunissen/particle_swarm
    cd particle_swarm
    make

### Author
Jannis Teunissen, jannis@teunissen.net
