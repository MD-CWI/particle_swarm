Electron swarm simulator
====

This code can be used to perform electron swarm simulations, using a Monte Carlo
particle model. The Monte Carlo particle model is written in Fortran 90, and
uses some Fortran 2003 features. There is a Python 2/3 command line interface
which can be used to perform parallel swarm computations, see the example below.

How the transport data is computed is described in e.g. http://dx.doi.org/10.1088/0022-3727/38/16/032

Author: Jannis Teunissen, jannis@teunissen.net

### Requirements

1. Particle model -- gfortran 4.8 or newer
2. Command line interface -- Python 2.7 or newer with numpy

The code has been tested on GNU/Linux and Mac.

### Getting the code

    git clone https://github.com/jannisteunissen/particle_swarm

### Compiling the code

    cd particle_swarm
    make

### Using the command line interface

    # To see command line options
    ./swarm_cli.py --help

    # Example of invocation
    ./swarm_cli.py input/cs_example.txt results.txt -gc N2 1.0 -flin 1e7 2e7 10

The last command uses cross sections from the file input/cs_example.txt, stores
results in results.txt and performs swarm simulations in nitrogen (N2) for 10
electric fields between 1e7 V/m and 2e7 V/m.

### Without the command line interface

    # Compute transport data in a single electric field,
    # with the settings from config.txt
    ./particle_swarm config.txt

A file config_example.txt is included, which shows the basic options. In the
output directory (specified by output_dir in the configuration file), a file
[...]_config.txt will appear, which shows all available options with some
documentation.

### Getting input data (cross sections)

Input data can be obtained from http://lxcat.net, the downloaded text files with
cross sections can directly be used.

### Format of output data

The first line of the output is a header with the names of the columns. Below
the header there is a matrix of values, separated by spaces. Each row gives data
for one electric field.

### Related software / projects

* Bolos: https://github.com/aluque/bolos
* Bolsig+: http://www.bolsig.laplace.univ-tlse.fr
* LXCat: http://lxcat.net
* Magboltz: http://consult.cern.ch/writeup/magboltz/
