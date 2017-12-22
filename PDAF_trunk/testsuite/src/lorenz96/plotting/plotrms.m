function[] = plotrms(varargin)
% plot_rms('path and file' [,further_arguments])
%
% Opens a NetCDF holding assimilation output from assimilating 
% into the Lorenz96 model and plots the true and estimated
% rms errors. 
%
% 'path and file': File name including path
%
% 'further_arguments': 'plot_forecast': 1 for yes (default)
%                      'plot_analysis': 1 for yes (default)
%
% This file is part of the test suite of PDAF.

global estimated_rmse_ini 
global true_rmse_ini
global mean_estimated_rmse_for
global mean_estimated_rmse_ana
global mean_true_rmse_for
global mean_true_rmse_ana

global trmse_for
global trmse_ana
global rmse_for
global rmse_ana

if length(varargin)==0
  disp('Function arguments incomplete - see help!')
  return
end

filename = varargin{1}

plot_forecast = 1;
plot_analysis = 1;
plot_meanerror = 0;
if length(varargin)>1
   for i=2:length(varargin):2
      arg = varargin{i};
      argval = varargin{i+1};
      if (length(arg)==13)
         if arg=='plot_forecast'
            plot_forecast = argval;
         end
         if arg=='plot_analysis'
            plot_analysis = argval;
         end
      end
   end
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

% Read time steps
varid = netcdf.inqVarID(nc,'time');
time = netcdf.getVar(nc,varid);
varid = netcdf.inqVarID(nc,'step');
step = netcdf.getVar(nc,varid);

% Read errors
varid = netcdf.inqVarID(nc,'rmse_ini');
rmse_ini = netcdf.getVar(nc,varid);
varid = netcdf.inqVarID(nc,'trmse_ini');
trmse_ini = netcdf.getVar(nc,varid);

varid = netcdf.inqVarID(nc,'rmse_for');
rmse_for = netcdf.getVar(nc,varid);
varid = netcdf.inqVarID(nc,'trmse_for');
trmse_for = netcdf.getVar(nc,varid);

varid = netcdf.inqVarID(nc,'rmse_ana');
rmse_ana = netcdf.getVar(nc,varid);
varid = netcdf.inqVarID(nc,'trmse_ana');
trmse_ana = netcdf.getVar(nc,varid);

varid = netcdf.inqVarID(nc,'mrmse_for_null');
mrmse_for = netcdf.getVar(nc,varid);
varid = netcdf.inqVarID(nc,'mrmse_ana_null');
mrmse_ana = netcdf.getVar(nc,varid);

varid = netcdf.inqVarID(nc,'mtrmse_for_null');
mtrmse_for = netcdf.getVar(nc,varid);
varid = netcdf.inqVarID(nc,'mtrmse_ana_null');
mtrmse_ana = netcdf.getVar(nc,varid);

netcdf.close(nc);

% Plot
if plot_analysis==1
   hf1 = figure;
   plot(rmse_ana,'k')
   hold
   plot(trmse_ana,'b')
   title(['Analysis errors for time step ',num2str(step(1)), ' to ', num2str(step(length(step)))]);
   legend('estimated error','true error')
end

if plot_forecast==1
   hf2 = figure;
   plot(rmse_for,'k')
   hold
   plot(trmse_for,'b')
   title(['Forecast errors for time step ',num2str(step(1)), ' to ', num2str(step(length(step)))]);
   legend('estimated error','true error')
end

% Display true and estimated mean errors
estimated_rmse_ini = rmse_ini 
true_rmse_ini = trmse_ini
mean_estimated_rmse_for = mrmse_for
mean_estimated_rmse_ana = mrmse_ana
mean_true_rmse_for = mtrmse_for
mean_true_rmse_ana = mtrmse_ana

