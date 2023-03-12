function[] = plot_obs(varargin)
% plot_obs('filename with path', timestep)
%
% Opens a NetCDF holding observations for the Lorenz63 model
% and plots the observation at a selected time.
%
% Arguments:
% 'filename with path': File name including path
% iteration           : Time step in file to be shown
%
% This file is part of the test suite of PDAF.

global obs;

if length(varargin)<2
  disp('Function arguments incomplete - see help!')
  return
end

% Name of file holding state trajectory
filename = varargin{1}

% Iteration in file to be shown
iter = varargin{2}-1

% Open file
if exist(filename,'file')
  nc=netcdf.open(filename,'nowrite');
  varid = netcdf.inqDimID(nc,'timesteps');
  [varname, n_steps] = netcdf.inqDim(nc, varid);

  disp(['file contains ',int2str(n_steps), ' timesteps'])
else
  disp('file does not exist!')
end

% Read state dimension
varid = netcdf.inqDimID(nc,'dim_state');
[varname dim] = netcdf.inqDim(nc,varid);

% Read state
varid = netcdf.inqVarID(nc,'obs');
obs = netcdf.getVar(nc,varid,[0,iter],[dim,1]);
varid = netcdf.inqVarID(nc,'time');
time = netcdf.getVar(nc,varid,iter,1);
varid = netcdf.inqVarID(nc,'step');
step = netcdf.getVar(nc,varid,iter,1);

netcdf.close(nc);

% Plot state
hf=figure;
plot(obs,'b+-')
title(['Observations for Lorenz63 model at time ',num2str(time),' (time step ',num2str(step),')'])
