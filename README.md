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

Some of the important options:

* Select particle mover: `-mover {analytic,boris,verlet}`
* Electric field range (in V/m): `-E_range min max`
* Number of electric fields: `-E_num number`
* How to vary electric fields: `-E_vary {lin,log}`
* Similar options for the magnetic field (in Tesla): `-B_range`, `-B_num`, `-B_vary`
* Similar options for angle between E and B (in degrees): `-angle_range`, `-angle_num`, `-angle_vary`

Note that by default `E_num = 1`, `B_num = 1`, and `angle_num = 1`.

## More examples

Compute transport coefficients between `E = 1e6 V/m` and `E = 2e7 V/m`, using a
logarithmic spacing with 20 points:

    ./swarm_cli.py input/cs_example.txt -of test.txt -gc N2 1.0 -E_range 1e6 2e7
        -E_num 20 -E_vary log

Compute transport coefficients in an electric field between `1e7 V/m` and `2e7
V/m`, for a magnetic field of `20 T` and `E,B`-angles between 5 and 85 degrees.
Use linearly spaced electric field and angles (the default), with 10 points:

    ./swarm_cli.py input/cs_example.txt -of test.txt -gc N2 1.0 -E_range 1e7 2e7
        -E_num 10 -B_range 20 20 -angle_range 5 85 -angle_num 10


Compute transport coefficients in an electric field between `1e7 V/m` and `2e7
V/m`, for a perpendicular magnetic field between `5 T` and `20 T`, using
linearly spaced electric and magnetic fields (the default) with 10 points:

    ./swarm_cli.py input/cs_example.txt -of test.txt -gc N2 1.0 -E_range 1e7 2e7
        -E_num 10 -angle_range 90 90 -B_range 5 20 -B_num 10

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

1. `Bz`: magnetic field (Tesla)
2. `Ez`: z-component of electric field (V/m)
3. `Ey`: y-component of electric field (V/m)
4. `angle`: angle between **E** and **B** in degrees
5. `omega_c`: gyration frequency (rad/s)
6. `energy`: mean energy (eV)
7. `mu_E`: mobility parallel to **E** (m2/Vs)
8. `mu_B`: mobility parallel to **B** (m2/Vs)
9. `mu_xB`: mobility perpendicular to **B** (m2/Vs)
10. `mu_ExB`: mobility in **ExB**-direction (m2/Vs) (the **ExB**-velocity over `Ey`)
11. `alpha`: ionization coefficient (1/m)
12. `eta`: attachment coefficient (1/m)
13. `coll_rate`: actual collision rate (1/s)
14. `diff_1`: x-component of diagonal diffusion coefficient (m2/s)
15. `diff_2`: y-component
16. `diff_3`: z-component
17. `vel_1`: x-component of mean velocity (m/s)
18. `vel_2`: y-component
19. `vel_3`: z-component
20. `vel_sq_1`: x-component of mean squared velocity
21. `vel_sq_2`: y-component
22. `vel_sq_3`: z-component

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
