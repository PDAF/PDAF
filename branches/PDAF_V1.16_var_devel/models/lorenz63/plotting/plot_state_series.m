function[] = plot_state_series(varargin)
% plot_state('filename with path', variable [, choice])
%
% Opens NetCDF output from the Lorenz63 model
% and plots sa selected state variable over time.
%
% Arguments:
% 'filename with path': File name including path
% variable            : variable to show
% choice              : Type of state to plot
%       choices: t - true, f - forecast, a - analysis
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

if length(varargin)>2
  plottype = varargin{3}
end

% Variable to show
vari = varargin{2}-1

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
time = netcdf.getVar(nc,varid,0,n_steps);
varid = netcdf.inqVarID(nc,'step');
step = netcdf.getVar(nc,varid,0,n_steps);

% Read state
if plottype=='t'
  varid = netcdf.inqVarID(nc,'state');
  state = netcdf.getVar(nc,varid,[vari,0],[1,n_steps]);
  statestr = 'true state';
else
  if plottype=='f'
    varid = netcdf.inqVarID(nc,'state_for');
    statestr = 'forecast estimate';
  elseif plottype=='a'
    varid = netcdf.inqVarID(nc,'state_ana');
    statestr = 'analysis estimate';
  end
  state = netcdf.getVar(nc,varid,[vari,0],[1,n_steps]);
end

netcdf.close(nc);

if vari==1
  strvari = 'x';
elseif vari==2
  strvari = 'y';
elseif vari==3
  strvari = 'z';
end

% Plot state
hf=figure;
plot(time, state,'r')
ylabel(strvari)
xlabel('time')  
title(['Lorenz63 model ' statestr ' variable ' strvari])



 
