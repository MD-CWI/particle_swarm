# Electron swarm simulator

This is a code for electron swarm simulations, using a Monte Carlo particle
model. Simply said, electrons are repeatedly released in a homogeneous and
static **E,B** field configuration. After they have relaxed to the background
field their properties are recorded, from which transport coefficient such as
the drift velocity or the ionization coefficient can be computed.

The key features of this code are:

* Support for an arbitrary **E,B** configuration
* Automatic variance/error estimation
* Automatic resizing of electrons swarms when there is ionization/attachment
* Support for different particle movers (Verlet, Boris pusher, analytic scheme)
* Written in modern Fortran
* A Python 2/3 command line interface for parallel computations

How the transport data can be computed from electron properties is described in
e.g. [Petrovic et al. 2005](http://dx.doi.org/10.1088/0022-3727/38/16/032).

People involved:

* Jannis Teunissen: original author
* Anbang Sun: testing, fixing bugs, first implementation of magnetic field

## Requirements

1. Particle model: gfortran 4.8 or newer
2. Command line interface: Python 2.7 or newer with numpy

The code has been tested on GNU/Linux and Mac.

## Getting the code

    git clone https://gitlab.com/MD-CWI-NL/particle_swarm

## Compiling the code

    cd particle_swarm
    make

## Using the command line interface

To see command the line options and their documentation, use:

    ./swarm_cli.py --help

An example of an invocation:

    ./swarm_cli.py input/cs_example.txt -of test.txt -gc N2 1.0 -E_range 1e7 2e7 -E_num 10

The last command uses cross sections from the file `input/cs_example.txt`,
stores results in `test.txt` and performs swarm simulations in pure nitrogen
(N2) for 10 electric fields. The fields are linearly interpolated between `1e7
V/m` and `2e7 V/m`.

A run can use a single electric field or a range of values, specified in `V/m`,
for example as:

    -E 1e7
    -E_range 1e7 2e7 -E_num 10

For a range of values, the spacing can be specified as `-E_vary lin` (linear,
the default) or `-E_vary log` (logarithmic). Similarly, a single magnetic field
or a range of magnetic fields can be specified in Tesla, using for example:

    -B 10
    -B_range 10 20 -B_num 10

By default, the magnetic field is zero. Finally, the angle (in degrees) between
the electric and magnetic can be varied:

    -angle 90
    -angle_range 10 80 -angle_num 5

By default, the electric and magnetic field are parallel. For the magnetic field
and, the spacing can be specified similarly as for the electric field,
with `-B_vary`.

Three different particle movers have been implemented, which can be selected
using `-mover {analytic,boris,verlet}`.

## More examples

Compute transport coefficients between `E = 1e6 V/m` and `E = 2e7 V/m`, using a
logarithmic spacing with 20 points:

    ./swarm_cli.py input/cs_example.txt -gc N2 1.0 -E_range 1e6 2e7 -E_num 20 -E_vary log

Compute transport coefficients in an electric field between `1e7 V/m` and `2e7
V/m`, for a magnetic field of `20 T` and `E,B`-angles between 5 and 85 degrees.
Use linearly spaced electric field and angles (the default), with 10 points:

    ./swarm_cli.py input/cs_example.txt -gc N2 1.0 -E_range 1e7 2e7 -E_num 10 -B 20 -angle_range 5 85 -angle_num 10

Compute transport coefficients in an electric field between `1e7 V/m` and `2e7
V/m`, for a perpendicular magnetic field between `5 T` and `20 T`, using
linearly spaced electric and magnetic fields (the default) with 10 points:

    ./swarm_cli.py input/cs_example.txt -gc N2 1.0 -E_range 1e7 2e7 -E_num 10 -angle 90 -B_range 5 20 -B_num 10

## Converting tables

The resulting table with transport data can be converted to a format that is
perhaps more useful in a discharge simulation with the conversion scripts in the
`tools` directory.

## Usage without the command line interface

    # Compute transport data with the settings from config.txt
    ./particle_swarm config.txt

A file `config_example.txt` is included, which shows some of the basic options.
In the output directory (specified by `output_dir` in the configuration file), a
file `[...]_config.txt` will appear, which shows all available options with some
documentation.

## Getting input data (cross sections)

Input data can be obtained from http://lxcat.net, the downloaded text files with
cross sections can directly be used.

## Format of output data

The first line of the output is a header with the names of the columns, below
which there is a matrix of values. This matrix contains the following columns:

1. `E`: electric field (V/m)
2. `B`: magnetic field (Tesla)
3. `angle`: angle between **E** and **B** in degrees. `B` always points in the z-direction, and `E` lies in the y,z-plane, so that `Ez = cos(angle) * E` and `Ey = sin(angle) * E`
4. `omega_c`: gyration frequency (rad/s)
5. `energy`: mean energy (eV)
6. `mu_E`: mobility parallel to **E** (m2/Vs)
7. `mu_B`: mobility parallel to **B** (m2/Vs)
8. `mu_xB`: mobility perpendicular to **B** (m2/Vs)
9. `mu_ExB`: mobility in **ExB**-direction (m2/Vs) (the **ExB**-velocity over `Ey`)
10. `alpha`: ionization coefficient (1/m)
11. `eta`: attachment coefficient (1/m)
12. `coll_rate`: actual collision rate (1/s)
13. `diff_1`: x-component of diagonal diffusion coefficient (m2/s)
14. `diff_2`: y-component
15. `diff_3`: z-component
16. `vel_1`: x-component of mean velocity (m/s)
17. `vel_2`: y-component
18. `vel_3`: z-component
19. `vel_sq_1`: x-component of mean squared velocity
20. `vel_sq_2`: y-component
21. `vel_sq_3`: z-component

## Which particle mover to use

Currently, three particle movers are included, which can be selected with the
`-mover <name>` option.

* `verlet`: The simplest mover, which can be used when there is no magnetic
  field.
* `boris`: Use Boris' method to push particles. Can be used for arbitrary
  electric and magnetic fields. The method slows down for very high magnetic
  fields because it requires time steps shorter than the gyration time.
* `analytic`: Update the particle position and velocity analytically.

## A comment on performance

How long a swarm needs to run depends on two factors:

1. The required accuracy: as with most Monte Carlo methods, an `N` times longer
   run gives a `sqrt(N)` times better accuracy.
2. The relaxation time of electrons. In general, electrons relax slower to the
   background field at lower energies.

## TODO

* Perform comparison with other solvers
* Include measurement of 'bulk' transport data
* Think of smart way to use individual particle trajectories

## Related software / projects

* Bolos: https://github.com/aluque/bolos
* Bolsig+: http://www.bolsig.laplace.univ-tlse.fr
* LXCat: http://lxcat.net
* Magboltz: http://consult.cern.ch/writeup/magboltz/
* METHES: http://fr.lxcat.net/download/METHES/
