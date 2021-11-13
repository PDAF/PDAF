function[] = plot_obs(varargin)
% plot_obs('filename with path', variable)
%
% Opens a NetCDF holding observations for the Lorenz63 model
% and plots the observation at a selected time.
%
% Arguments:
% 'filename with path': File name including path
% variable            : variable to show (1=x, 2=y, 3=z)
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
vari = varargin{2}-1

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
obs = netcdf.getVar(nc,varid,[vari,0],[1,n_steps]);
varid = netcdf.inqVarID(nc,'time');
time = netcdf.getVar(nc,varid,0,n_steps);
varid = netcdf.inqVarID(nc,'step');
step = netcdf.getVar(nc,varid,0,n_steps);

netcdf.close(nc);

% Set string for state variable
if vari==0
  strvari = 'x';
elseif vari==1
  strvari = 'y';
elseif vari==2
  strvari = 'z';
end

% Plot state
hf=figure;
plot(obs,'b')
title(['Observations for Lorenz63 model, variable ' strvari])
