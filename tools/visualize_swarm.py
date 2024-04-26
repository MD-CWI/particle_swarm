#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from glob import glob

electron_volt = 1.602176634e-19
electron_mass = 9.10938189e-31

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='''Visualize electron swarm''')
parser.add_argument('basename', type=str, default='visualize/part_data',
                    help='Base filename of output files')
parser.add_argument('-figsize', type=float, nargs=2, default=[8., 5.],
                    help='Figure size')
parser.add_argument('-save', type=str,
                    help='Save animation to file')
parser.add_argument('-i_x', type=int, default=2,
                    help='Spatial index to use on x-axis')
parser.add_argument('-i_y', type=int, default=0,
                    help='Spatial index to use on y-axis')
parser.add_argument('-max_eV', type=float, default=10.,
                    help='Maximum energy in eV for visualization')
args = parser.parse_args()

files = sorted(glob(args.basename + '*.txt'))

N = len(files)
positions = []
velocities = []
energies_eV = []
times = np.zeros(N)
x_min = np.ones(3) * 1e10
x_max = np.ones(3) * -1e10

for i, fname in enumerate(files):
    with open(fname) as f:
        first_line = f.readline().strip()
        times[i] = float(first_line.split()[1])
        coords = np.loadtxt(f)

        if coords.ndim == 1:
            coords = np.reshape(coords, [1, -1])

        r = coords[:, 0:3] * 1e3
        v = coords[:, 3:6]
        eV = 0.5 * electron_mass * np.sum(v**2, axis=1) / electron_volt

        x_min = np.minimum(x_min, r.min(axis=0))
        x_max = np.maximum(x_max, r.max(axis=0))

        positions.append(r.T)
        velocities.append(v.T)
        energies_eV.append(eV.T)


plt.rcParams['font.size'] = 12

fig = plt.figure(figsize=(args.figsize), layout='tight')
dx = 0.1 * (x_max[args.i_x] - x_min[args.i_x])
dy = 0.1 * (x_max[args.i_y] - x_min[args.i_y])
ax = fig.add_subplot(autoscale_on=False,
                     xlim=(x_min[args.i_x] - dx, x_max[args.i_x] + dx),
                     ylim=(x_min[args.i_y] - dx, x_max[args.i_y] + dx))
ax.set_aspect('equal')
ax.set_xlabel('mm')
ax.set_ylabel('mm')

scat = ax.scatter([], [], c=[], alpha=1.0, s=20., vmin=0.,
                  vmax=args.max_eV, marker='o', cmap='cool')
time_text = ax.text(0.1, 0.8, '', transform=ax.transAxes)
plt.colorbar(scat, ax=ax, location='top', label='eV', shrink=0.5)


def animate(i):
    scat.set_offsets(np.c_[positions[i][args.i_x],
                           positions[i][args.i_y]])
    scat.set_color(scat.to_rgba(energies_eV[i]))

    n_particles = len(positions[i][args.i_x])
    time_text.set_text(f't = {times[i]*1e9:.2f} ns, N = {n_particles}')
    return scat, time_text


ani = animation.FuncAnimation(
    fig, animate, N, interval=50, blit=True)

if args.save:
    ani.save(args.save, fps=30, dpi=200)
else:
    plt.show()
