#!/usr/bin/env python

# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import tempfile
import shutil
import os
import sys
import numpy as np
from subprocess import check_output, CalledProcessError


def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Command line interface to compute electron
        transport data from swarm simulations.
        Author: Jannis Teunissen, jannis@teunissen.net''',
        epilog='''Usage example: ./swarm_cli.py input/cs_example.txt
        -gc N2 1.0 -E_range 1e7 2e7 -E_num 10''')
    parser.add_argument('cs', type=str,
                        help='File with cross sections')
    parser.add_argument('-name', type=str, default='visualize/part_data',
                        help='Name (including path) of output files')
    parser.add_argument('-video', action='store_true',
                        help='Try to generate a video (requires gnuplot + ffpmeg)')
    parser.add_argument('-plot',
                        help='Custom gnuplot command')
    parser.add_argument('-fps', type=float, default=10.,
                        help='Frames per second for animation')
    parser.add_argument('-gc', dest='gas_comps', type=str, nargs='+',
                        required=True, metavar='gas frac',
                        help='List of gas names and fractions, '
                        'for example: N2 0.8 O2 0.2')
    parser.add_argument('-mover', type=str,
                        choices=['analytic', 'boris', 'verlet'],
                        default='analytic', help='Choice of particle mover')

    parser.add_argument('-t', type=float, required=True,
                        help='End time for visualization')
    parser.add_argument('-dt', type=float, required=True,
                        help='Time step for visualization')
    parser.add_argument('-npart', type=int, default=10,
                        help='Initial number of particles')
    parser.add_argument('-maxpart', type=int, default=1000000,
                        help='Maximum number of particles')
    parser.add_argument('-E', type=float, required=True,
                        help='Electric field value (V/m)')
    parser.add_argument('-angle', type=float, default=0.0,
                        help='angle between E and B (degrees)')
    parser.add_argument('-B', type=float, default=0.0,
                        help='Magnetic field (T)')

    parser.add_argument('-T', type=float, default=300.,
                        metavar='temperature', help='Gas temperature (K)')
    parser.add_argument('-p', type=float, default=1.0,
                        metavar='pressure', help='Gas pressure (bar)')
    return parser.parse_args()


def create_swarm_cfg(tmpdir, args):
    fname = tmpdir + '/base_cfg.txt'
    f = open(fname, 'w')
    f.write('output_dir = ' + tmpdir + '\n')
    f.write('gas_file = ' + args.cs + '\n')
    f.write('swarm_name = visualize\n')
    f.write('gas_components = ' + ' '.join(args.gas_comps[0::2]) + '\n')
    f.write('gas_fractions = ' + ' '.join(args.gas_comps[1::2]) + '\n')
    f.write('gas_pressure = ' + str(args.p) + '\n')
    f.write('gas_temperature = ' + str(args.T) + '\n')
    f.write('swarm_size = ' + str(args.npart) + '\n')
    f.write('consecutive_run = F\n')
    f.write('dry_run = F\n')
    f.write('visualize_only = T\n')
    f.write('visualize_rotate_Ez = T\n')
    f.write('visualize_end_time = ' + str(args.t) + '\n')
    f.write('visualize_dt_output = ' + str(args.dt) + '\n')
    f.write('visualize_max_particles = ' + str(args.maxpart) + '\n')
    f.write('visualize_base_name = ' + args.name + '\n')
    f.write('electric_field = ' + str(args.E) + '\n')
    f.write('magnetic_field = ' + str(args.B) + '\n')
    f.write('field_angle_degrees = ' + str(args.angle) + '\n')
    f.write('particle_mover = ' + args.mover + '\n')
    f.close()
    return fname


def pswarm_wrapper(cmd_and_num):
    cmd, num = cmd_and_num
    try:
        res = check_output(cmd)
    except:
        print("Error when running: " + cmd)
        sys.exit(1)
    num.value += 1
    return res


def print_swarm_info(args):
    print("Starting particle swarm visualization")
    print("----------------------------------------")
    print("Cross section file : {}".format(args.cs))
    print("Gas components     : {}".format(' '.join(args.gas_comps[0::2])))
    print("Gas fractions      : {}".format(' '.join(args.gas_comps[1::2])))
    print("Temperature (K)    : {}".format(args.T))
    print("Pressure (bar)     : {}".format(args.p))
    print("Particle mover     : {}".format(args.mover))
    print("Output base name   : {}".format(args.name))
    print("----------------------------------------")
    print("Electric field (V/m) : {}".format(args.E))
    print("Magnetic field (T)   : {}".format(args.B))
    print("Angle (degrees)      : {}".format(args.angle))
    print("----------------------------------------")


if __name__ == '__main__':
    args = get_args()

    print_swarm_info(args)

    try:
        tmpdir = tempfile.mkdtemp(dir=os.getcwd())

        args.dirname, args.basename = os.path.split(args.name)
        if args.dirname == '':
            args.dirname = '.'

        if not os.path.exists(args.dirname):
            os.makedirs(args.dirname)

        base_cfg = create_swarm_cfg(tmpdir, args)

        # Perform a dry run of particle_swarm to generate lookup tables for the
        # cross sections
        try:
            out = check_output(['./particle_swarm', base_cfg])
        except CalledProcessError as e:
            print("particle_swarm returned an error")
            print(e.output)
            sys.exit(1)

        gpscript = args.dirname + '/script.gp'
        ix_max = int(round(args.t/args.dt))
        last_file = '{}_{:06d}.txt'.format(args.name, ix_max)

        with open(last_file) as f:
            lines = f.readlines()
            # Sensitive to file syntax!
            rmin = np.array(map(float, lines[2].split()[6:9]))
            rmax = np.array(map(float, lines[3].split()[6:9]))

        dsize = rmax - rmin
        rmax = rmin + dsize.max()

        if args.video:
            with open(gpscript, 'w') as gp:
                gp.write("set terminal pngcairo size 800,600 enhanced\n" +
                         "unset xtics; set xlabel 'x'\n" +
                         "unset ytics; set ylabel 'y'\n" +
                         "unset ztics; set zlabel 'z'\n" +
                         "set xrange [{}:{}]\n".format(rmin[0], rmax[0]) +
                         "set yrange [{}:{}]\n".format(rmin[1], rmax[1]) +
                         "set zrange [{}:{}]\n".format(rmin[2], rmax[2]) +
                         "set view equal xyz\n" +
                         "do for [i=0:{}] {{\n".format(ix_max) +
                         "set output sprintf('{}_%06d.png', i)\n".format(args.basename) +
                         "splot sprintf('{}_%06d.txt', i)".format(args.basename) +
                         " u 1:2:3 notitle w p pt 7 ps 1}\n")

            os.system('cd {} && gnuplot script.gp'.format(args.dirname))
            cmd = 'ffmpeg -y -r {} -f image2 -i {}_%06d.png -vframes {} {}.mp4'.format(
                args.fps, args.name, ix_max+1, args.name)
            os.system(cmd)
            print("Generated {}.mp4".format(args.name))

    finally:
        shutil.rmtree(tmpdir)

