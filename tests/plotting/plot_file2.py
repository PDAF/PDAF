#!/usr/bin/env python3

# A script to plot the 2D field from the online_2D_serialmodel tutorial.
# Requires Python 3, Matplotlib and Numpy.

# Usage: ./plot_field.py <filename>

import matplotlib.pyplot as plt
import numpy as np
import argparse as ap

def read_and_plot(filename):
    field = np.loadtxt(filename)
    field = field.reshape(512,2048)
    plt.imshow(field, origin='lower',interpolation='none')
    plt.title(filename,fontsize=8)
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
