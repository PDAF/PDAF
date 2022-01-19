#!/usr/bin/env python3

# A script to plot a time series of the ensemble mean and spread
# for a single point from the online_2D_serialmodel tutorial.
# Requires Python 3, Matplotlib and Numpy.

# Usage: ./plot_ens.py <i> <j>


import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import glob
import sys

def read_ensemble():
    files = glob.glob("ens_01*.txt")
    ntime = len(files)
    files = glob.glob("ens_*step02_for.txt")
    nens = len(files)

    try:
        if nens*ntime == 0:
            raise Exception("No ens*.txt files found.")
    except:
        raise

    ens = np.empty((ntime, 18, 36, nens))
    for i in range(int(ntime/2)):
        for j in range(nens):
            field = np.loadtxt("ens_{0:02d}_step{1:02d}_for.txt".format(j+1, 2*(i+1)))
            ens[2*i, ..., j] = np.loadtxt("ens_{0:02d}_step{1:02d}_for.txt".format(j+1, 2*(i+1)))
            ens[2*i+1, ..., j] = np.loadtxt("ens_{0:02d}_step{1:02d}_ana.txt".format(j+1, 2*(i+1)))
    return ens

def plot(ens, i, j):
    time = np.array([])
    for k in range(int(ens.shape[0]/2)):
        time = np.hstack((time, [k, k]))
    plt.fill_between(time, ens[:, i, j].mean(axis=-1)-ens[:, i, j].std(axis=-1), ens[:, i, j].mean(axis=-1)+ens[:, i, j].std(axis=-1), color=[0.7,0.7,0.7])
    plt.plot(time, ens[:, i, j].mean(axis=-1))
    plt.title('Ensemble at '+str(i)+ ', '+str(j))
    plt.gca().set_xticklabels(["{0:.0f}".format(i+2) for i in 2*plt.gca().get_xticks()])
    plt.show()

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument('i', type=int)
    parser.add_argument('j', type=int)
    args = parser.parse_args()
    print("Plotting time series for point at i={0:d}, j={0:d}.".format(args.i, args.j))
    try:
        ens = read_ensemble()
    except Exception as err:
        print(err)
        sys.exit(1)

    plot(ens, args.i, args.j)
