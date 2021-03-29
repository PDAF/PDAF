function[] = plot_state_series(varargin)
% plot_state('file1', 'file2', variable [, choice1, choice2])
%
% Opens NetCDF output from the Lorenz63 model
% and plots two time series (e.g. truth and analysis).
%
% Arguments:
% 'file1': Name of file 1 including path
% 'file2': Name of file 2 including path
% variable            : variable to show (1:x, 2:y, 3:z)
% choice1              : Type of state to plot (default: t)
% choice2              : Type of state to plot (default: a)
%       choices: t - true, f - forecast, a - analysis
%
% This file is part of the test suite of PDAF.

% Default is to plot the true state
plottype1 = 't';
plottype2 = 'a';

if length(varargin)<3
  disp('Function arguments incomplete - see help!')
  return
end

% Name of file holding state trajectory
filename1 = varargin{1}
filename2 = varargin{2}
vari = varargin{3}-1

if length(varargin)>3
  plottype1 = varargin{4}
end
if length(varargin)>4
  plottype2 = varargin{5}
end

% Iteration in file to be shown
if plottype1~='t'
    vari1 = varargin{3}-1
else
    vari1 = varargin{3}
end
if plottype2~='t'
    vari2 = varargin{3}-1
else
    vari2 = varargin{3}
end

% Open files
%if exist(filename1,'file1')
  nc1=netcdf.open(filename1,'nowrite');
  varid1 = netcdf.inqUnlimDims(nc1);
  [varname1, n_steps1] = netcdf.inqDim(nc1, varid1);

  disp(['file 1 contains ',int2str(n_steps1), ' timesteps'])    
%else
%  disp('file 1 does not exist!')
%end
%if exist(filename2,'file2')
  nc2 = netcdf.open(filename2,'nowrite');
  varid2 = netcdf.inqUnlimDims(nc2);
  [varname2, n_steps2] = netcdf.inqDim(nc2, varid2);

  disp(['file 2 contains ',int2str(n_steps2), ' timesteps'])    
%else
%  disp('file 2 does not exist!')
%end

% Read state dimension
varid = netcdf.inqDimID(nc1,'dim_state');
[varname1 dim] = netcdf.inqDim(nc1,varid);

% Read time and time step
varid = netcdf.inqVarID(nc1,'time');
time1 = netcdf.getVar(nc1,varid,0,n_steps1);
varid = netcdf.inqVarID(nc1,'step');
step1 = netcdf.getVar(nc1,varid,0,n_steps1);
varid = netcdf.inqVarID(nc2,'time');
time2 = netcdf.getVar(nc2,varid,0,n_steps2);
varid = netcdf.inqVarID(nc2,'step');
step2 = netcdf.getVar(nc2,varid,0,n_steps2);

% Read states
if plottype1=='t'
  varid = netcdf.inqVarID(nc1,'state');
  state1 = netcdf.getVar(nc1,varid,[vari,0],[1,n_steps1]);
  statestr1 = 'true state';
else
  if plottype1=='f'
    varid = netcdf.inqVarID(nc1,'state_for');
    statestr1 = 'forecast estimate';
  elseif plottype1=='a'
    varid = netcdf.inqVarID(nc1,'state_ana');
    statestr1 = 'analysis estimate';
  end
  state1 = netcdf.getVar(nc1,varid,[vari,0],[1,n_steps1]);
end
if plottype2=='t'
  varid = netcdf.inqVarID(nc2,'state');
  state2 = netcdf.getVar(nc2,varid,[vari,0],[1,n_steps2]);
  statestr2 = 'true state';
else
  if plottype2=='f'
    varid = netcdf.inqVarID(nc2,'state_for');
    statestr2 = 'forecast estimate';
  elseif plottype2=='a'
    varid = netcdf.inqVarID(nc2,'state_ana');
    statestr2 = 'analysis estimate';
  end
  state2 = netcdf.getVar(nc2,varid,[vari,0],[1,n_steps2]);
end
netcdf.close(nc1);
netcdf.close(nc2);

if vari==1
  strvari = 'x';
elseif vari==2
  strvari = 'y';
elseif vari==3
  strvari = 'z';
end

% Plot states
hf=figure;
plot(time1, state1,'b')
hold on
plot(time2, state2,'k')
ylabel(strvari)
xlabel('time') 
legend(statestr1, statestr2)
title(['Lorenz63 model, variable ' strvari])



 
