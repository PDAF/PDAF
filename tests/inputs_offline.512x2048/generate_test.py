# Script to generate input files for the PDAF offline tutorial
# This is a translation of the former Matlab script to Python
#
# Note, that the random numbers are different in Matlab and
# Python, so that the varification outputs need to be consistent.
#     
# L. Nerger, 9/2024

import numpy as np
      
dim_x = 2048        # Grid dimension in x-direction
dim_y = 512        # Grid dimension in y-direction
dim_ens = 9        # Maximum ensemble size
stddev_obs = 0.5   # Observation error standard deviation
dxobs = 10          # x-Grid spacing for observations type A
dyobs = 9          # y-Grid spacing for observations type A
dxobsB = 17        # x-Grid spacing for observations type B
dyobsB = 16         # y-Grid spacing for observations type B
obsB_offsetx = -4  # x-offset in position of observations type B
obsB_offsety = -2  # y-offset in position of observations type B
#voidx_min=12
#voidx_max=25
      
dowrite = 1        # 1 to write files


# Initialize random number generator
np.random.seed(0)

# Generate true field

field = np.zeros((dim_y, dim_x))
for j in range(dim_x):
   for i in range(dim_y):
      field[i,j] = np.sin(np.pi*((i+1)/dim_y + (j+1)/dim_x))

if dowrite==1:
   np.savetxt('true.txt', field[:,:])


# Generate ensemble states

ens = np.zeros((dim_y, dim_x, dim_ens))
for k in range(dim_ens):
   for j in range(dim_x):
      for i in range(dim_y):
         ens[i,j,k] = np.sin(np.pi*((i+1)/dim_y + (j+1)/dim_x) + 0.5*np.pi*(k+1+5)/dim_ens)

if dowrite==1:
   for k in range(dim_ens):
      np.savetxt('ens_'+str(k+1)+'.txt', ens[:,:,k])


# Compute ensemble mean = forecast state

state = np.mean(ens,axis=2)
if dowrite==1:
   np.savetxt('state_ini.txt', state[:,:])


# Observations

obs_error = np.zeros((dim_y, dim_x))
full_obs = np.zeros((dim_y, dim_x))
obs_error = stddev_obs * np.random.randn(dim_y, dim_x) 

full_obs[:,:] = field[:,:] + obs_error


# Obs type A
obs = np.zeros((dim_y, dim_x))
obs[:,:] = -999

for j in range(dxobs-1,dim_x,dxobs):
   for i in range(dyobs-1,dim_y,dyobs):
      obs[i,j] = full_obs[i,j]

# Introduce a gap in x-direction to have dim_obs_p=0 with parallelization
#obs[:,voidx_min:voidx_max] = -999

if dowrite==1:
   np.savetxt('obs.txt', obs[:,:])

# Obs type B
obs = np.zeros((dim_y, dim_x))
obs[:,:] = -999

for j in range(dxobsB-1+obsB_offsetx,dim_x,dxobsB):
   for i in range(dyobsB-1+obsB_offsety,dim_y,dyobsB):
      obs[i,j] = full_obs[i,j]

if dowrite==1:
   np.savetxt('obsB.txt', obs[:,:])

