#!/usr/bin/env python

"""
Plot functions for the Lorenz96 example of the PDAF test
suite. This file is part of the test suite of PDAF.

To use as a Python script, pass the function name followed by the
parameters as command line arguments. Also see ./plot_l96.py --help.

08/2018 Gernot Geppert, University of Reading

"""


import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import argparse as ap


def plot_example():
    """
    Plot time-mean RMS errors for the set of example experiments with
    the Lorenz96 model provided in tools/runasml.sh.

    Plotted is the time-mean true rms analysis error for assimilation
    with the ESTKF filter over 10000 time steps. The observations are
    used from time step 1000 onwards. We vary the forgetting
    factor. Here, the filter diverges for forgetting factors above
    0.98.

    """

    forgets = np.arange(1, 0.89, -0.01)
    mean_rmse = np.zeros_like(forgets)
    mean_trmse = np.zeros_like(forgets)

    for i in range(len(forgets)):
        with nc.Dataset('../t1_N30_f{0:g}.nc'.format(forgets[i])
                        ) as ncfile:
            mrmse_for  = ncfile['mrmse_for_null'][0]
            mtrmse_for = ncfile['mtrmse_for_null'][0]
            mrmse_ana  = ncfile['mrmse_ana_null'][0]
            mtrmse_ana = ncfile['mtrmse_ana_null'][0]

        mean_rmse[i] = mrmse_ana
        mean_trmse[i] = mtrmse_ana

    fig, ax = plt.subplots()

    ax.plot(forgets, mean_trmse, 'b+-', label="true error")
    ax.plot(forgets, mean_rmse, 'k+-', label="estimated error")

    ax.legend()
    ax.set_xlabel("forgetting factor")
    ax.set_ylabel("mean RMS error")
    ax.set_title("Time-mean true RMS analysis errors for ESTKF, N=30")


def plot_eofs(filename, index):
    """
    Open a NetCDF file holding an EOF-decomposed covariance matrix for the
    Lorenz96 model and plots the selected eigenvectors  or the mean state.

    Parameters
    ----------
    filename : str
        file name including path
    index : int
        index of eigenvector to plot (0 for mean state)

    """

    with nc.Dataset(filename) as ncfile:
        print("File contains {0:d} eigenvectors.".format(
                ncfile.dimensions['rank'].size
                                                         )
              )
        if index > 0:
            state = ncfile['u_svd'][index-1]
        else:
            state = ncfile['meanstate'][0]

    fig, ax = plt.subplots()
    ax.plot(range(1, len(state)+1), state, 'b')
    if index > 0:
        ax.set_title("Eigenvector {0:d} of Lorenz96 model".format(index))
    else:
        ax.set_title("Mean state of Lorenz96 model")


def plot_obs(filename, timestep):
    """
    Open a NetCDF file holding observations for the Lorenz96 model and
    plots the observation at a selected time.

    Parameters
    ----------
    filename : str
        file name including path
    timestep : int
        time step to plot

    """

    with nc.Dataset(filename) as ncfile:
        print("File contains {0:d} time steps.".format(
                ncfile.dimensions['timesteps'].size
                                                       )
              )
        obs  = ncfile['obs'][timestep-1]
        time = ncfile['time'][timestep-1]
        step = ncfile['step'][timestep-1]

    fig, ax = plt.subplots()
    ax.plot(range(1, len(obs)+1), obs, 'b+-')
    ax.set_title(("Observations for Lorenz96 model at time {0:.3f} "
                  "(time step {1:d})").format(time, step)
                 )

def plot_rms(filename, plot_forecast=False, plot_analysis=True):
    """
    Open a NetCDF holding assimilation output from assimilating into
    the Lorenz96 model and plos the true and estimated rms errors.

    Parameters
    ----------
    filename : str
        file name including path
    plot_forecast : bool, optional
        plot forecast RMSEs
    plot_analysis : bool, optional
        plot analysis RMSEs

    """

    with nc.Dataset(filename) as ncfile:
        print("File contains {0:d} time steps.".format(
                ncfile.dimensions['iteration'].size
                                                       )
              )
        step      = ncfile['step'][:]
        rmse_ini  = ncfile['rmse_ini'][0]
        trmse_ini = ncfile['trmse_ini'][0]
        rmse_for  = ncfile['rmse_for'][:]
        trmse_for = ncfile['trmse_for'][:]
        rmse_ana  = ncfile['rmse_ana'][:]
        trmse_ana = ncfile['trmse_ana'][:]
        mrmse_for  = ncfile['mrmse_for_null'][0]
        mtrmse_for = ncfile['mtrmse_for_null'][0]
        mrmse_ana  = ncfile['mrmse_ana_null'][0]
        mtrmse_ana = ncfile['mtrmse_ana_null'][0]

    if plot_forecast:
        fig_forecast, ax_forecast = plt.subplots()
        ax_forecast.plot(range(1, len(rmse_for)+1), rmse_for, 'k', label="estimated error")
        ax_forecast.plot(range(1, len(trmse_for)+1), trmse_for, 'b', label="true error")
        ax_forecast.set_title(("Forecast errors for time step {0:d} to {1:d}"
                               "").format(step[0], step[-1])
                              )
        ax_forecast.legend()

    if plot_analysis:
        fig_analysis, ax_analysis = plt.subplots()
        ax_analysis.plot(range(1, len(rmse_ana)+1), rmse_ana, 'k', label="estimated error")
        ax_analysis.plot(range(1, len(trmse_ana)+1), trmse_ana, 'b', label="true error")
        ax_analysis.set_title(("Analysis errors for time step {0:d} to {1:d}"
                               "").format(step[0], step[-1])
                              )
        ax_analysis.legend()

    print("Estimated initial RMSE:      ", np.mean(rmse_ini))
    print("True initial RMSE:           ", np.mean(trmse_ini))
    print("Estimated mean forecast RMSE:", np.mean(mrmse_for))
    print("True mean forecast RMSE:    ", np.mean(mtrmse_for))
    print("Estimated mean analysis RMSE:", np.mean(mrmse_ana))
    print("True mean analysis RMSE:     ", np.mean(mtrmse_ana))


def plot_sigma(filename):
    """
    Open a NetCDF file holding an EOF-decomposed covariance matrix for
    the Lorenz96 model and plot the eigenvalue spectrum.

    Parameters
    ----------
    filename : str
        file name including path
    timestep : int
        time step to plot

    """

    with nc.Dataset(filename) as ncfile:
        sigma = ncfile['sigma'][:]

    fig, ax = plt.subplots()
    ax.plot(range(1, len(sigma)+1), sigma, 'b')
    ax.set_title("Eigenvalues of covariance matrix for Lorenz96 model")


def plot_state(filename, timestep, choice='t'):
    """
    Open NetCDF output from the Lorenz96 model and plot the state at
    the selected time.

    Parameters
    ----------
    filename : str
        file name including path
    timestep : int
        time step to plot
    choice : str
        't' for true state, 'f' for forecast, 'a' for analysis,
        'i' for initial state. Default is 't'.

    """

    if choice not in 'tfai':
        raise Exception(("{0} is not a valid choice. Valid choices are "
                         "'t', 'f', 'a', 'i'.".format(choice))
                        )

    if choice == 'i':
        timestep = 0
    elif choice == 't':
        timestep = timestep
    else:
        timestep = timestep-1

    with nc.Dataset(filename) as ncfile:
        print("File contains {0:d} time steps.".format(
                ncfile['step'].size
                                                       )
              )
        time = ncfile['time'][timestep]
        step = ncfile['step'][timestep]

        if choice == 't':
            state    = ncfile['state'][timestep]
            statestr = 'true state'
        elif choice == 'i':
            state    = ncfile['state_ini'][timestep]
            statestr = 'initial state'
        elif choice == 'f':
            state    = ncfile['state_for'][timestep]
            statestr = 'forecast estimate'
        else:
            state    = ncfile['state_ana'][timestep]
            statestr = 'analysis estimate'

    fig, ax = plt.subplots()
    ax.plot(range(1, len(state)+1), state, 'r')
    ax.set_title(("Lorenz96 model {0} at time {1:.3f} (time step {2:d})"
                  "").format(statestr, time, step)
                 )


def plot_diff(file1, file2, timestep, choice1='t', choice2='a'):
    """
    Open NetCDF output from the Lorenz96 model and plot the difference
    of the state at the selected time.

    Parameters
    ----------
    file1 : str
        file name 1 including path
    file2 : str
        file name 2 including path
    timestep : int
        time step to plot
    choice1 : str
    choice2 : str
        't' for true state, 'f' for forecast, 'a' for analysis,
        'i' for initial state. Default is 't'.

    """

    if choice1 not in 'tfai':
        raise Exception(("{0} is not a valid choice. Valid choices are "
                         "'t', 'f', 'a', 'i'.".format(choice))
                        )
    if choice2 not in 'tfai':
        raise Exception(("{0} is not a valid choice. Valid choices are "
                         "'t', 'f', 'a', 'i'.".format(choice))
                        )

    if choice1 == 'i':
        timestep1 = 0
    elif choice1 == 't':
        timestep1 = timestep
    else:
        timestep1 = timestep-1

    if choice2 == 'i':
        timestep2 = 0
    elif choice1 == 't':
        timestep2 = timestep
    else:
        timestep2 = timestep-1

    with nc.Dataset(file1) as ncfile1:
        print("File contains {0:d} time steps.".format(
                ncfile1['step'].size
                                                       )
              )
        time1 = ncfile1['time'][timestep1]
        step1 = ncfile1['step'][timestep1]

        if choice1 == 't':
            state1    = ncfile1['state'][timestep1]
            statestr1 = 'true state'
        elif choice1 == 'i':
            state1    = ncfile1['state_ini'][timestep1]
            statestr1 = 'initial state'
        elif choice1 == 'f':
            state1    = ncfile1['state_for'][timestep1]
            statestr1 = 'forecast estimate'
        else:
            state1    = ncfile1['state_ana'][timestep1]
            statestr1 = 'analysis estimate'

    with nc.Dataset(file2) as ncfile2:
        print("File contains {0:d} time steps.".format(
                ncfile2['step'].size
                                                       )
              )
        time2 = ncfile2['time'][timestep2]
        step2 = ncfile2['step'][timestep2]

        if choice2 == 't':
            state2    = ncfile2['state'][timestep2]
            statestr2 = 'true state'
        elif choice2 == 'i':
            state2    = ncfile2['state_ini'][timestep2]
            statestr2 = 'initial state'
        elif choice2 == 'f':
            state2    = ncfile2['state_for'][timestep2]
            statestr2 = 'forecast estimate'
        else:
            state2    = ncfile2['state_ana'][timestep2]
            statestr2 = 'analysis estimate'

    diff = state1 - state2
    
    fig, ax = plt.subplots()
    ax.plot(range(1, len(state1)+1), diff, 'r')
    ax.set_title(("Difference at time {1:.3f} (time step {2:d})"
                  "").format(statestr1, time1, step2)
                 )


def _parse_str(arg):
    """
    Convert string to integer or False, if possible.

    'False' is converted to the bool constant False, 'True' does not
    need to be converted as every non-empty string evaluates as True.

    Parameters
    ----------
    arg : str
        string to be checked and possibly converted

    Returns
    -------
    str, int or bool

    """

    try:
        return int(arg)
    except ValueError:
        if arg == 'False':
            return False
        else:
            return arg

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument('function', type=str,
                        help="name of the function to be used"
                        )
    parser.add_argument('arguments', nargs='*', metavar='argument',
                        help=("arguments passed to function "
                              "(see plot_l96.py for documentation)")
                        )

    args = parser.parse_args()
    args.arguments = [_parse_str(arg) for arg in args.arguments]

    locals()[args.function](*args.arguments)

    plt.show()
