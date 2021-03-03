function[] = plotsigma(varargin)
% plotsigma('filename with path')
%
% Opens a NetCDF holding an EOF-decomposed covariance matrix 
% for the Lorenz96 model and plots the eigenvalue spectrum
%
% Arguments:
% 'filename with path': File name including path
%
% This file is part of the test suite of PDAF.

if length(varargin)<1
  disp('Function arguments incomplete - see help!')
  return
end

% Name of file holding state trajectory
filename = varargin{1}

% Open File
if exist(filename,'file')
  nc=netcdf.open(filename,'nowrite');
else
  disp('file does not exist!')
end

% Read singular values
varid = netcdf.inqVarID(nc,'sigma');
sigma = netcdf.getVar(nc,varid);

netcdf.close(nc);

% Plot singular values
hf=figure;
plot(sigma,'b')
title(['Eigenvalues of covariance matrix for Lorenz96 model'])




 
