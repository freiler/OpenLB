#!/usr/bin/env python3

"""
Post-processing sample, written in Python

Copyright (C) 2020 Louis Kronberg, Stephan Simonis
E-mail contact: info@openlb.net
The most recent release of OpenLB can be downloaded at
<http://www.openlb.net/>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public
License along with this program; if not, write to the Free
Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
Boston, MA  02110-1301, USA.

./ade2dplot.py

This script creates plots from csv files generated by the ade2d.cpp simulation program.
More specifically this script:
    - plots a loglog plot with reference curves is created from
      ./tmp/gnuplotData/data/averageSimL2RelErr.dat
    - plots the error between analytical solution and numerical solution for each timestep

The script expects the csv files to be organized in a specific Directory structure.
If the name of the output files/directories are changed then this script likely has to
be changed aswell.

Note: that changing the possible values of the discretisation paramter N
does not require any changes to this script. This script assumes all csv files to be in ./tmp.
This behaviour can be changed by changing the value of the OUTPUT_PATH global variable.
"""

import sys
import csv
from pathlib import Path
from math import nan, inf
import matplotlib.pyplot as plt

if sys.version_info[0] < 3:
    print("ERROR: " + __file__ + " must be run with python3")
    sys.exit(1)

# specify all files that should be plot in a loglog plot
LOGLOG_PATH = Path('./tmp/gnuplotData/data/averageSimL2RelErr.dat')
OUTPUT_PATH = Path('./tmp/')

# proportionality constants for the reference curves in the loglog plot
C1 = 10
C2 = 100

extract_int = lambda string: int(''.join(s for s in string if s.isdigit()))

def read_datafile(filename):
    x = []
    y = []
    with open(filename, 'r') as datfile:
        data = csv.reader(datfile, delimiter=' ')
        for row in data:
            x.append(float(row[0]))
            y.append(float(row[-1]))
    if nan in x or nan in y:
        print("[Warning]: {filename} contains Nan's")
    if inf in x or inf in y:
        print("[Warning]: {filename} contains Inf's")
    return x, y

def plot(filepath, **kwargs):
    outfilename = kwargs.pop('outfilename', filepath.resolve().parents[1].joinpath(filepath.stem + '.png'))
    saveplot = kwargs.pop('saveplot', True)
    kwargs.setdefault('title', filepath.stem)
    kwargs.setdefault('label', filepath.stem)

    x, y = read_datafile(filepath)
    if not x and not y:
        return None

    plot_func = kwargs.pop('plot_func', plt.plot)
    ylim = kwargs.pop('ylim', None)
    if ylim:
        plt.ylim(*ylim)

    plt.title(kwargs.pop('title', ''))
    plt.xlabel(kwargs.pop('xlabel', ''))
    plt.ylabel(kwargs.pop('ylabel', ''))
    plt.margins(*kwargs.pop('margin', (0,0)))
    plot_func(x, y, **kwargs)
    plt.legend()
    plt.grid(kwargs.pop('grid', True), which='both')

    if saveplot:
        plt.savefig(outfilename)
        plt.close('all')


if __name__ == "__main__":
    simdirs = (path for path in OUTPUT_PATH.iterdir() if path.is_dir() and path.name.startswith('N'))

    if not OUTPUT_PATH.exists():
        sys.exit(f"Abort. The output path: {OUTPUT_PATH} does not exist. You have to run the simulation first before using this script")

    if LOGLOG_PATH.exists():
        x, y = read_datafile(LOGLOG_PATH)
        xref = range(int(min(x)), int(max(x)))
        plt.loglog(xref, [C1*x**-1 for x in xref], label='$\mathcal{O}(N^{-1})$')
        plt.loglog(xref, [C2*x**-2 for x in xref], label='$\mathcal{O}(N^{-2})$')
        plot(LOGLOG_PATH,
             plot_func=plt.loglog,
             margin=(0.1,0.1),
             marker='o',
             linestyle='',
             fillstyle='none',
             xlabel='N',
             title='Experimental order of convergence', label='time-averaged rel L2 err')
    else:
        print(f"Skipping log-log plot with reference curves. Because {LOGLOG_PATH} does not exist")


    # iterate over each simulation's tmp folder
    for simdir in simdirs:
        data_dir = simdir / 'gnuplotData' / 'data'
        for filename in data_dir.iterdir():
            if filename.suffix == '.dat':
                plot(filename, 
		     label='L2 err', 
		     margin=(.1,.1),
                     xlabel = 'timestep')
