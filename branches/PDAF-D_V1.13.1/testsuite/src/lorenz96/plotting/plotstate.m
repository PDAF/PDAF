function[] = plotstate(varargin)
% plotstate('filename with path', iteration [, type])
%
% Opens NetCDF output from the Lorenz96 model
% and plots teh state at a selected time.
%
% Arguments:
% 'filename with path': File name including path
% iteration           : Interation in file to be shown
% type                : Type of state to plot
%       choices: t - true, f - forecast, a - analysis, i - inial
%
% This file is part of the test suite of PDAF.

% Default is to plot the true state
plottype = 't';

if length(varargin)<2
  disp('Function arguments incomplete - see help!')
  return
end

% Name of file holding state trajectory
filename = varargin{1}

% Iteration in file to be shown
iter = varargin{2}-1

if length(varargin)>2
  plottype = varargin{3}
end

if plottype=='i'
  iter = 0
end

% Open file
if exist(filename,'file')
  nc=netcdf.open(filename,'nowrite');
  varid = netcdf.inqUnlimDims(nc);
  [varname, n_steps] = netcdf.inqDim(nc, varid);

  disp(['file contains ',int2str(n_steps), ' timesteps'])    
else
  disp('file does not exist!')
end

% Read state dimension
varid = netcdf.inqDimID(nc,'dim_state');
[varname dim] = netcdf.inqDim(nc,varid);

% Read time and time step
varid = netcdf.inqVarID(nc,'time');
time = netcdf.getVar(nc,varid,iter,1);
varid = netcdf.inqVarID(nc,'step');
step = netcdf.getVar(nc,varid,iter,1);

% Read state
if plottype=='t'
  varid = netcdf.inqVarID(nc,'state');
  state = netcdf.getVar(nc,varid,[0,iter],[dim,1]);
  statestr = 'true state';
else
  if plottype=='i'
    varid = netcdf.inqVarID(nc,'state_ini');
    statestr = 'initial state';
  elseif plottype=='f'
    varid = netcdf.inqVarID(nc,'state_for');
    statestr = 'forecast estimate';
  elseif plottype=='a'
    varid = netcdf.inqVarID(nc,'state_ana');
    statestr = 'analysis estimate';
  end
  state = netcdf.getVar(nc,varid,[0,iter],[dim,1]);
end

netcdf.close(nc);

% Plot state
hf=figure;
plot(state,'r')
if plottype=='i'
  title(['Lorenz96 model ' statestr ' at time ',num2str(time),' (time step ',num2str(step),')'])
else
  title(['Lorenz96 model ' statestr ' at time ',num2str(time),' (time step ',num2str(step),')'])
end



 
