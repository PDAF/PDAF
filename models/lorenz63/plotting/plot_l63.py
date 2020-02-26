#!/usr/bin/env python

"""
Plot functions for the Lorenz63 example of the PDAF test
suite. This file is part of the test suite of PDAF.

To use as a Python script, pass the function name followed by the
parameters as command line arguments. Also see ./plot_l63.py --help.

07/2019 Lars Nerger based on Script for Lorenz86 by Gernot Geppert, University of Reading

"""


import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import argparse as ap


def plot_eofs(filename, index):
    """
    Open a NetCDF file holding an EOF-decomposed covariance matrix for the
    Lorenz63 model and plots the selected eigenvectors  or the mean state.

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
        ax.set_title("Eigenvector {0:d} of Lorenz63 model".format(index))
    else:
        ax.set_title("Mean state of Lorenz63 model")


def plot_obs(filename, timestep):
    """
    Open a NetCDF file holding observations for the Lorenz63 model and
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
    ax.set_title(("Observations for Lorenz63 model at time {0:.3f} "
                  "(time step {1:d})").format(time, step)
                 )


def plot_obs_series(filename, variable):
    """
    Open a NetCDF file holding observations for the Lorenz63 model and
    plots a time series for the selected variable.

    Parameters
    ----------
    filename : str
        file name including path
    variable : int
        variable to plot (1=x, 2=y, 3=z)

    """

    with nc.Dataset(filename) as ncfile:
        print("File contains {0:d} time steps.".format(
                ncfile.dimensions['timesteps'].size
                                                       )
              )
        obs  = ncfile['obs'][:,variable-1]
        time = ncfile['time'][:]
        step = ncfile['step'][:]

    if variable==1:
        strvari = 'X'
    elif variable==2:
        strvari = 'Y'
    elif variable==3:
        strvari = 'Z'

    fig, ax = plt.subplots()
    ax.plot(time, obs, 'b')
    ax.set_title(("Observations for Lorenz63 model, variable {0}".format(strvari))
                 )

def plot_rms(filename, plot_forecast=False, plot_analysis=True):
    """
    Open a NetCDF holding assimilation output from assimilating into
    the Lorenz63 model and plos the true and estimated rms errors.

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
    the Lorenz63 model and plot the eigenvalue spectrum.

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
    ax.set_title("Eigenvalues of covariance matrix for Lorenz63 model")


def plot_state(filename, timestep, choice='t'):
    """
    Open NetCDF output from the Lorenz63 model and plot the state at
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
    ax.set_title(("Lorenz63 model {0} at time {1:.3f} (time step {2:d})"
                  "").format(statestr, time, step)
                 )



def plot_state_series(filename, variable, choice='t'):
    """
    Open NetCDF output from the Lorenz63 model and plot a time
    series of the selected variable.

    Parameters
    ----------
    filename : str
        file name including path
    variable : int
        variable to plot (1=x, 2=y, 3=z)
    choice : str
        't' for true state, 'f' for forecast, 'a' for analysis,
        'i' for initial state. Default is 't'.

    """

    if choice not in 'tfai':
        raise Exception(("{0} is not a valid choice. Valid choices are "
                         "'t', 'f', 'a', 'i'.".format(choice))
                        )

    with nc.Dataset(filename) as ncfile:
        print("File contains {0:d} time steps.".format(
                ncfile['step'].size
                                                       )
              )
        time = ncfile['time'][:]
        step = ncfile['step'][:]

        if choice == 't':
            state    = ncfile['state'][:,variable-1]
            statestr = 'true state'
        elif choice == 'i':
            state    = ncfile['state_ini'][:,variable-1]
            statestr = 'initial state'
        elif choice == 'f':
            state    = ncfile['state_for'][:,variable-1]
            statestr = 'forecast estimate'
        else:
            state    = ncfile['state_ana'][:,variable-1]
            statestr = 'analysis estimate'

    if variable==1:
        strvari = 'X'
    elif variable==2:
        strvari = 'Y'
    elif variable==3:
        strvari = 'Z'

    fig, ax = plt.subplots()
    ax.plot(time, state, 'r')
    ax.set_title(("Lorenz63 model {0}, variable {1}"
                  "").format(statestr, strvari)
                 )


def plot_2state_series(filename1, filename2, variable, choice1='t', choice2='a'):
    """
    Open NetCDF output from the Lorenz63 model and plot two time
    series of the selected variable.

    Parameters
    ----------
    filename1 : str
        file name including path
    filename2 : str
        file name including path
    variable : int
        variable to plot (1=x, 2=y, 3=z)
    choice1 : str  (choice for file 1)
        't' for true state, 'f' for forecast, 'a' for analysis,
        'i' for initial state. Default is 't'.
    choice2 : str  (choice for file 2)
        't' for true state, 'f' for forecast, 'a' for analysis,
        'i' for initial state. Default is 't'.

    """

    if choice1 not in 'tfai':
        raise Exception(("{0} is not a valid choice. Valid choices are "
                         "'t', 'f', 'a', 'i'.".format(choice1))
                        )
                        
    if choice2 not in 'tfai':
        raise Exception(("{0} is not a valid choice. Valid choices are "
                         "'t', 'f', 'a', 'i'.".format(choice2))
                        )

    with nc.Dataset(filename1) as ncfile1:
        print("File contains {0:d} time steps.".format(
                ncfile1['step'].size
                                                       )
              )
        time1 = ncfile1['time'][:]
        step1 = ncfile1['step'][:]

        if choice1 == 't':
            state1    = ncfile1['state'][:,variable-1]
            statestr1 = 'true state'
        elif choice1 == 'i':
            state1    = ncfile1['state_ini'][:,variable-1]
            statestr1 = 'initial state'
        elif choice1 == 'f':
            state1    = ncfile1['state_for'][:,variable-1]
            statestr1 = 'forecast estimate'
        else:
            state1    = ncfile1['state_ana'][:,variable-1]
            statestr1 = 'analysis estimate'

    with nc.Dataset(filename2) as ncfile2:
        print("File contains {0:d} time steps.".format(
                ncfile2['step'].size
                                                       )
              )
        time2 = ncfile2['time'][:]
        step2 = ncfile2['step'][:]

        if choice2 == 't':
            state2    = ncfile2['state'][:,variable-1]
            statestr2 = 'true state'
        elif choice2 == 'i':
            state2    = ncfile2['state_ini'][:,variable-1]
            statestr2 = 'initial state'
        elif choice1 == 'f':
            state2    = ncfile2['state_for'][:,variable-1]
            statestr2 = 'forecast estimate'
        else:
            state2    = ncfile2['state_ana'][:,variable-1]
            statestr2 = 'analysis estimate'

    if variable==1:
        strvari = 'X'
    elif variable==2:
        strvari = 'Y'
    elif variable==3:
        strvari = 'Z'

    fig, ax = plt.subplots()
    ax.plot(time1, state1, 'r',label=statestr1)
    ax.plot(time2, state2, 'b',label=statestr2)
    ax.legend()
    ax.set_title(("Lorenz63 model, variable {0}"
                  "").format(strvari)
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
                              "(see plot_l63.py for documentation)")
                        )

    args = parser.parse_args()
    args.arguments = [_parse_str(arg) for arg in args.arguments]

    locals()[args.function](*args.arguments)

    plt.show()
