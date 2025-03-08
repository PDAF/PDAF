# Script to generate input files for the PDAF offline tutorial
# This is a translation of the former Matlab script to Python
#
# Note, that the random numbers are different in Matlab and
# Python, so that the varification outputs need to be consistent.
#     
# L. Nerger, 9/2024

import numpy as np
      
dim_x = 36         # Grid dimension in x-direction
dim_y = 18         # Grid dimension in y-direction
dim_ens = 9        # Maximum ensemble size
stddev_obs = 0.5   # Observation error standard deviation
dxobs = 5          # x-Grid spacing for observations type A
dyobs = 4          # y-Grid spacing for observations type A
dxobsB = 11        # x-Grid spacing for observations type B
dyobsB = 8         # y-Grid spacing for observations type B
obsB_offsetx = -4  # x-offset in position of observations type B
obsB_offsety = -2  # y-offset in position of observations type B
      
dowrite = 1        # 1 to write files

# Locations of observations not placed at grid points (x, y)
obs_interp = np.zeros((11,2))
obs_interp[:,:] = [[3.0, 2.1],
     [3.4, 6.8],
     [6.1, 6.8],
     [8.9, 7.6],
     [8.9, 14.9],
     [20.0, 6.4],
     [20.4, 16.1],
     [14.1, 10.2],
     [31.0, 5.2],
     [31.2, 11.9],
     [28.9, 14.9]];

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


# Compute ensemble mean = initial state estimate

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


# Interpolated observations

iobs_error = np.zeros((dim_y, dim_x))
iobs_error = stddev_obs * np.random.randn(len(obs_interp)) 

gx = np.zeros(2)
gy = np.zeros(2)
iobs = np.zeros((len(obs_interp),3))
for i in range(len(obs_interp)):
   # Get closest grid points
   gx[0] = np.floor(obs_interp[i,0])
   gx[1] = np.ceil(obs_interp[i,0])
   if gx[1]==gx[0]:
      gx[1] = gx[1]+1
   gy[0] = np.floor(obs_interp[i,1])
   gy[1] = np.ceil(obs_interp[i,1])
   if gy[1]==gy[0]:
      gy[1] = gy[1]+1

   # Compute interpolation coefficients
   icoeff = np.zeros(4)
   denum = (gx[1]-gx[0])*(gy[1]-gy[0]);
   icoeff[0] = (gx[1] - obs_interp[i,0]) * (gy[1] - obs_interp[i,1])/denum
   icoeff[1] = (obs_interp[i,0] - gx[0]) * (gy[1] - obs_interp[i,1])/denum
   icoeff[2] = (gx[1] - obs_interp[i,0]) * (obs_interp[i,1] - gy[0])/denum
   icoeff[3] = (obs_interp[i,0] - gx[0]) * (obs_interp[i,1] - gy[0])/denum

   # Interpolate
   iobs[i,0] = icoeff[0]*field[np.int(gy[0]-1),np.int(gx[0]-1)] + icoeff[1]*field[np.int(gy[0]-1),np.int(gx[1]-1)] + icoeff[2]*field[np.int(gy[1]-1),np.int(gx[0]-1)] + icoeff[3]*field[np.int(gy[1]-1),np.int(gx[1]-1)]

   # Add error
   iobs[i,0] = iobs[i,0] + iobs_error[i]

   # Augment with coordinates
   iobs[i,1] = obs_interp[i,0]
   iobs[i,2] = obs_interp[i,1]

if dowrite==1:
   obsfile = open(r"obsC.txt", "w")
   obsfile.write(str(len(obs_interp))+'\n')
   for i in range(len(obs_interp)):
      obsfile.write(str(iobs[i,0])+' '+str(iobs[i,1])+' '+str(iobs[i,2])+'\n')
   obsfile.close()


# Prepare full field of interpolated obs for plotting

if dowrite==1:
   iobs_full = np.zeros((dim_y,dim_x))
   iobs_full[:,:] = -999
   for i in range(len(obs_interp)):
      iobs_full[np.int(np.floor(obs_interp[i,1])), np.int(np.floor(obs_interp[i,0]))] = iobs[i,0]

   np.savetxt('obsC_field.txt', iobs_full[:,:])

