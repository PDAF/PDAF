function[] = ploteofs(varargin)
% ploteofs('filename with path', mode index)
%
% Opens a NetCDF holding an EOF-decomposed covariance matrix 
% for the Lorenz96 model and plots the selected eigenvectors,
% or the mean state.
%
% Arguments:
% 'filename with path': File name including path
% mode index          : Index of eigenvector to be shows (0 for mean state)
%
% This file is part of the test suite of PDAF.

if length(varargin)<2
  disp('Function arguments incomplete - see help!')
  return
end

% Name of file holding state trajectory
filename = varargin{1}

% Iteration in file to be shown
iter = varargin{2}

% Open file
if exist(filename,'file')
  nc=netcdf.open(filename,'nowrite');
  varid = netcdf.inqDimID(nc,'rank');
  [varname, n_modes] = netcdf.inqDim(nc, varid);

  disp(['file contains ',int2str(n_modes), ' eigenvectors'])
else
  disp('file does not exist!')
end

% Read state dimension
varid = netcdf.inqDimID(nc,'dim_state');
[varname dim] = netcdf.inqDim(nc,varid);

% Read state or singular vector
if iter>0
  varid = netcdf.inqVarID(nc,'u_svd');
  state = netcdf.getVar(nc,varid,[0,iter-1],[dim,1]);
else
  varid = netcdf.inqVarID(nc,'meanstate');
  state = netcdf.getVar(nc,varid,[0,0],[dim,1]);
end

netcdf.close(nc);

% Plot state
hf=figure;
plot(state,'b')
if iter>0
   title(['Eigenvector ', num2str(iter), ' of Lorenz96 model'])
else
   title(['Meanstate of Lorenz96 model'])
end




 
