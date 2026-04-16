#!/usr/bin/env python3

# A script to plot the 2D field of observations from the online_2D_serialmodel tutorial.
# Requires Python 3, Matplotlib and Numpy.
# The script is identical to plot_field.py except for a color range limit.

# Usage: ./plot_obs.py <filename>

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import argparse as ap

def read_and_plot(filename):
    field = np.loadtxt(filename)
    field = field.reshape(18,36)
    field_m = ma.masked_outside(field,-10.0,10.0)
    plt.imshow(field_m, origin='lower',interpolation='none')
    plt.title(filename)
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()
    try:
        read_and_plot(args.filename)
    except OSError as err:
        print(err)
