#!/usr/bin/env python3

# A script to plot the difference of two 2D fields from the tutorial.
# Requires Python 3, Matplotlib and Numpy.
# This variant is for grid size 4: 512x512 grid points

# Usage: ./plot_diff4.py <filename1> <filename2>

import matplotlib.pyplot as plt
import numpy as np
import argparse as ap

def read_and_plot(filename1, filename2):
    field1 = np.loadtxt(filename1)
    field2 = np.loadtxt(filename2)
    field1 = field1.reshape(512,512)
    field2 = field2.reshape(512,512)
    rmse = 0;
    for i in range(512):
       for j in range(512):
	       rmse = rmse + (field1[i,j]-field2[i,j])**2
    rmse = np.sqrt(1/(512*512)*(rmse));
    print('RMSE: ', rmse)

    title = (filename1+' - '+filename2)
    plt.imshow(field1-field2, origin='lower',interpolation='none')
    plt.colorbar()
    plt.title(title,fontsize=8)
    plt.show()

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument('filename1')
    parser.add_argument('filename2')
    args = parser.parse_args()
    try:
        read_and_plot(args.filename1, args.filename2)
    except OSError as err:
        print(err)
