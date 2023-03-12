% plot_example
%
% Plot time-mean RMS errors for the set of example experiments
% with the Lorenz96 model provided in tools/runasml.sh.
%
% Plotted is the time-mean true rms analysis error for 
% assimilation with the SEIK filter over 10000 time steps. The
% observations are used from time step 1000 onwards. We vary the 
% forgetting factor. Here, the filter diverges for forgetting
% factors above 0.98.
%
% This file is part of the test suite of PDAF.

global mtrmse_ana;
global mrmse_ana;

forgets =[1 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.9]

for f=1:length(forgets)
   filename = ['../../../bin/t1_N30_f' num2str(forgets(f)) '.nc'];


   % Read errors
   nc=netcdf.open(filename,'nowrite');
   varid = netcdf.inqVarID(nc,'mrmse_for_null');
   mrmse_for = netcdf.getVar(nc,varid);
   varid = netcdf.inqVarID(nc,'mtrmse_for_null');
   mtrmse_for = netcdf.getVar(nc,varid);
   varid = netcdf.inqVarID(nc,'mrmse_ana_null');
   mrmse_ana = netcdf.getVar(nc,varid);
   varid = netcdf.inqVarID(nc,'mtrmse_ana_null');
   mtrmse_ana = netcdf.getVar(nc,varid);
   netcdf.close(nc);

   mean_trmse(f) = mtrmse_ana;
   mean_rmse(f) = mrmse_ana;
end

figure
plot(forgets, mean_trmse,'b+-');
hold
plot(forgets, mean_rmse,'k+-');
legend('true error','estimated error',2)
title('Time-mean true RMS analysis errors for SEIK, N=30')
xlabel('forgetting factor');
ylabel('mean RMS error');

