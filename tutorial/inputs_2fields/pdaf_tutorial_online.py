# Script to generate input files for the PDAF online tutorial for 2 fields;
# 'field' is the same field as in the one-fields case; 'fieldB' is added here
#
# This is a translation of the former Matlab script to Python
#
# Note, that the random numbers are different in Matlab and
# Python, so that the varification outputs need to be consistent.
#
# L. Nerger, 9/2024
# update for netCDF files, L. Nerger, 2/2026

import numpy as np
import netCDF4 as nc

dim_x = 36         # Grid dimension in x-direction
dim_y = 18         # Grid dimension in y-direction
dim_ens = 20       # Maximum ensemble size
dim_step = 50      # Number of time steps
stddev_obs = 0.5   # error standard deviation for observations type A
stddev_obsB = 0.25 # error standard deviation for observations type B
dxobs = 5          # x-Grid spacing for observations type A
dyobs = 4          # y-Grid spacing for observations type A
obs_offsetx = -4   # x-offset in position of observations type A
obs_offsety = -3   # y-offset in position of observations type A
dxobsB = 6         # x-Grid spacing for observations type B
dyobsB = 5         # y-Grid spacing for observations type B
obsB_offsetx = -2  # x-offset in position of observations type B
obsB_offsety = -1  # y-offset in position of observations type B
rotate = 1         # 1 to rotate ensemble states, 0 to shift them

dowrite = 1        # 1 to write files
write_nc = 0       # 1 to write netCDF files

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

# Generate true field A

field = np.zeros((dim_y, dim_x, dim_step+1))
for j in range(dim_x):
   for i in range(dim_y):
      field[i,j,0] = np.sin(2*np.pi*((i+1)/dim_y + (j+1)/dim_x))

for step in range(1,dim_step+1):
   for i in range(dim_y-1):
      field[i+1,:,step] = field[i,:,step-1]
   field[0,:,step] = field[-1,:,step-1]

# Write truth files for fieldA
if dowrite==1:
   np.savetxt('trueA_initial.txt', field[:,:,0])

   for step in range(1,dim_step+1):
      if step<10:
         stepstr = '0'+str(step)
      else:
         stepstr = str(step)

      np.savetxt('trueA_step'+stepstr+'.txt', field[:,:,step])

   if write_nc==1:
      for step in range(0, dim_step+1):
         if step==0:
            ncfile = nc.Dataset('trueA_ini.nc',mode='w')
         else:
            if step<10:
               stepstr = '0'+str(step)
            else:
               stepstr = str(step)

               ncfile = nc.Dataset('trueA_step'+stepstr+'.nc',mode='w')
         
               xdim = ncfile.createDimension('dim_x', dim_x)
               ydim = ncfile.createDimension('dim_y', dim_y)
               timedim = ncfile.createDimension('step', 1)
               trueA = ncfile.createVariable('field',np.float64, ('step', 'dim_x', 'dim_y',))

               trueA[0, :,:] = np.transpose(field[:,:,step])
               ncfile.close()

# Generate ensemble states A and their mean

ens = np.zeros((dim_y, dim_x, dim_ens))
for k in range(dim_ens):
   for j in range(dim_x):
      for i in range(dim_y):
         if rotate==1:
            # rotating ensemble states
            ens[i,j,k] = np.sin(2*np.pi*((i+1)/dim_y + (j+1)/dim_x*(0.2*(k-dim_ens/2))))
         else:
            # shifting ensemble states
            ens[i,j,k] = np.sin(2*np.pi*((i+1)/dim_y + (j+1)/dim_x + 0.25*k*np.pi/dim_y))

# Write ensemble files for fieldA
if dowrite==1:

   for k in range(dim_ens):
      if k<9:
         ensstr = '0'+str(k+1)
      else:
         ensstr = str(k+1)
      np.savetxt('ensA_'+ensstr+'.txt', ens[:,:,k])

   if write_nc==1:
      for k in range(dim_ens):
         if k<9:
            ensstr = '0'+str(k+1)
         else:
            ensstr = str(k+1)
         
            ncfile = nc.Dataset('ensA_'+ensstr+'.nc',mode='w')
         
            xdim = ncfile.createDimension('dim_x', dim_x)
            ydim = ncfile.createDimension('dim_y', dim_y)
            timedim = ncfile.createDimension('step', 1)
            trueB = ncfile.createVariable('field',np.float64, ('step', 'dim_x', 'dim_y',))

            trueB[0, :,:] = np.transpose(ens[:,:,k])
      ncfile.close()

# Compute ensemble mean = initial state estimate

state = np.mean(ens,axis=2)

# Write ensemble mean for fieldA
if dowrite==1:
   np.savetxt('stateA_ini.txt', state[:,:])

   if write_nc==1:
      ncfile = nc.Dataset('stateA_ini.nc',mode='w')
         
      xdim = ncfile.createDimension('dim_x', dim_x)
      ydim = ncfile.createDimension('dim_y', dim_y)
      timedim = ncfile.createDimension('step', 1)
      trueB = ncfile.createVariable('field',np.float64, ('step', 'dim_x', 'dim_y',))

      trueB[0, :,:] = np.transpose(state[:,:])
      ncfile.close()

# Generate true field B

fieldB = np.zeros((dim_y, dim_x, dim_step+1))
for j in range(dim_x):
   for i in range(dim_y):
      fieldB[i,j,0] = np.sin(4*np.pi*((i+1)/dim_y - (j+1)/dim_x))

for step in range(1,dim_step+1):
   for i in range(dim_y-1):
      fieldB[i+1,:,step] = fieldB[i,:,step-1]
   fieldB[0,:,step] = fieldB[-1,:,step-1]

# Write truth files for fieldB
if dowrite==1:

   np.savetxt('trueB_initial.txt', fieldB[:,:,0])

   for step in range(1,dim_step+1):
      if step<10:
         stepstr = '0'+str(step)
      else:
         stepstr = str(step)

      np.savetxt('trueB_step'+stepstr+'.txt', fieldB[:,:,step])

   if write_nc==1:
      for step in range(0, dim_step+1):
         if step==0:
            ncfile = nc.Dataset('trueB_ini.nc',mode='w')
         else:
            if step<10:
               stepstr = '0'+str(step)
            else:
               stepstr = str(step)

               ncfile = nc.Dataset('trueB_step'+stepstr+'.nc',mode='w')
         
            xdim = ncfile.createDimension('dim_x', dim_x)
            ydim = ncfile.createDimension('dim_y', dim_y)
            timedim = ncfile.createDimension('step', 1)
            trueB = ncfile.createVariable('field',np.float64, ('step', 'dim_x', 'dim_y',))

            trueB[0, :,:] = np.transpose(fieldB[:,:,step])
      ncfile.close()


# Generate ensemble states B and their mean

ensB = np.zeros((dim_y, dim_x, dim_ens))
for k in range(dim_ens):
   for j in range(dim_x):
      for i in range(dim_y):
         ensB[i,j,k] = np.sin(4*np.pi*((i+1)/dim_y - (j+1)/dim_x) + 3*0.5*np.pi*(k+1+5)/dim_ens) + 0.75*np.cos(3*np.pi*((i+1)/dim_y - (j+1)/dim_x))

# Write ensemble files for fieldB
if dowrite==1:

   for k in range(dim_ens):
      if k<9:
         ensstr = '0'+str(k+1)
      else:
         ensstr = str(k+1)
      np.savetxt('ensB_'+ensstr+'.txt', ensB[:,:,k])

   if write_nc==1:
      for k in range(dim_ens):
         if k<9:
            ensstr = '0'+str(k+1)
         else:
            ensstr = str(k+1)

         ncfile = nc.Dataset('ensB_'+ensstr+'.nc',mode='w')
         
         xdim = ncfile.createDimension('dim_x', dim_x)
         ydim = ncfile.createDimension('dim_y', dim_y)
         timedim = ncfile.createDimension('step', 1)
         trueB = ncfile.createVariable('field',np.float64, ('step', 'dim_x', 'dim_y',))

         trueB[0, :,:] = np.transpose(ensB[:,:,k])
      ncfile.close()

# Compute ensemble mean = initial state estimate

stateB = np.mean(ensB,axis=2)


# Write ensemble mean for fieldA
if dowrite==1:
   np.savetxt('stateB_ini.txt', stateB[:,:])

   if write_nc==1:
      ncfile = nc.Dataset('stateB_ini.nc',mode='w')
         
      xdim = ncfile.createDimension('dim_x', dim_x)
      ydim = ncfile.createDimension('dim_y', dim_y)
      timedim = ncfile.createDimension('step', 1)
      trueB = ncfile.createVariable('field',np.float64, ('step', 'dim_x', 'dim_y',))

      trueB[0, :,:] = np.transpose(stateB[:,:])
      ncfile.close()



# Observations A

obs_error = np.zeros((dim_y, dim_x, dim_step+1))
full_obs = np.zeros((dim_y, dim_x, dim_step+1))
obs_error = stddev_obs * np.random.randn(dim_y, dim_x, dim_step+1) 

full_obs[:,:,:] = field[:,:,:] + obs_error

obs = np.zeros((dim_y, dim_x, dim_step+1)) 
obs[:,:,:] = -999

for step in range(1,dim_step+1):
   for j in range(dxobs-1+obs_offsetx,dim_x,dxobs):
      for i in range(dyobs-1+obs_offsety,dim_y,dyobs):
         obs[i,j,step] = full_obs[i,j,step]

if dowrite==1:
# Option to write each observation time in a separate file
#    for step in range(1,dim_step+1):
#       ncfile = nc.Dataset('obsA_step'+str(step)+'.nc',mode='w')
#       xdim = ncfile.createDimension('dim_x', dim_x)
#       ydim = ncfile.createDimension('dim_y', dim_y)
#       timedim = ncfile.createDimension('step', 1)
#       true = ncfile.createVariable('obs',np.float64, ('step', 'dim_x', 'dim_y',))

#       true[0,:,:] = np.transpose(obs[:,:,step])
#    ncfile.close()

   for step in range(1,dim_step+1):
      if step<10:
         stepstr = '0'+str(step)
      else:
         stepstr = str(step)
      np.savetxt('obsA_step'+stepstr+'.txt', obs[:,:,step])

   if write_nc==1:
      # Write all observations into one file
      ncfile = nc.Dataset('obsA.nc',mode='w')
      xdim = ncfile.createDimension('dim_x', dim_x)
      ydim = ncfile.createDimension('dim_y', dim_y)
      timedim = ncfile.createDimension('step', dim_step)
      true = ncfile.createVariable('obs',np.float64, ('step', 'dim_x', 'dim_y',))

      for step in range(1,dim_step+1):
         true[step-1,:,:] = np.transpose(obs[:,:,step])
      ncfile.close()
         
# Observations B

obs_errorB = np.zeros((dim_y, dim_x, dim_step+1))
full_obsB = np.zeros((dim_y, dim_x, dim_step+1))
obs_errorB = stddev_obsB * np.random.randn(dim_y, dim_x, dim_step+1) 

full_obsB[:,:,:] = fieldB[:,:,:] + obs_errorB


obsB = np.zeros((dim_y, dim_x, dim_step+1)) 
obsB[:,:,:] = -999

for step in range(1,dim_step+1):
   for j in range(dxobsB-1+obsB_offsetx,dim_x,dxobsB):
      for i in range(dyobsB-1+obsB_offsety,dim_y,dyobsB): 
         obsB[i,j,step] = full_obsB[i,j,step]

if dowrite==1:
   for step in range(1,dim_step+1):
      if step<10:
         stepstr = '0'+str(step)
      else:
         stepstr = str(step)
      np.savetxt('obsB_step'+stepstr+'.txt', obsB[:,:,step])

   if write_nc==1:
      # Write all observations into one file
      ncfile = nc.Dataset('obsB.nc',mode='w')
      xdim = ncfile.createDimension('dim_x', dim_x)
      ydim = ncfile.createDimension('dim_y', dim_y)
      timedim = ncfile.createDimension('step', dim_step)
      true = ncfile.createVariable('obs',np.float64, ('step', 'dim_x', 'dim_y',))

      for step in range(1,dim_step+1):
         true[step-1,:,:] = np.transpose(obsB[:,:,step])
      ncfile.close()




# Interpolated observations

iobs_error = np.zeros((dim_y, dim_x, dim_step+1))
iobs_error = stddev_obs * np.random.randn(len(obs_interp), dim_step+1) 

gx = np.zeros(2)
gy = np.zeros(2)
iobs = np.zeros((len(obs_interp),3,dim_step+1))
for step in range(1,dim_step+1):
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
      iobs[i,0,step] = icoeff[0]*field[int(gy[0]-1),int(gx[0]-1),step] + icoeff[1]*field[int(gy[0]-1),int(gx[1]-1),step] + icoeff[2]*field[int(gy[1]-1),int(gx[0]-1),step] + icoeff[3]*field[int(gy[1]-1),int(gx[1]-1),step]

      # Add error
      iobs[i,0,step] = iobs[i,0,step] + iobs_error[i, step]

      # Augment with coordinates
      iobs[i,1,step] = obs_interp[i,0]
      iobs[i,2,step] = obs_interp[i,1]

if dowrite==1:
   for step in range(1,dim_step+1):
      obsfile = open(r"obsAint_step"+str(step)+".txt", "w")
      obsfile.write(str(len(obs_interp))+'\n')
      for i in range(len(obs_interp)):
         obsfile.write(str(iobs[i,0,step])+' '+str(iobs[i,1,step])+' '+str(iobs[i,2,step])+'\n')
      obsfile.close()

# Prepare full field for plotting

if dowrite==1:
   iobs_full = np.zeros((dim_y,dim_x,dim_step+1))
   iobs_full[:,:,:] = -999
   for step in range(1,dim_step+1):
      for i in range(len(obs_interp)):
         iobs_full[int(np.floor(obs_interp[i,1])), int(np.floor(obs_interp[i,0])),step] = iobs[i,0,step]

   for step in range(1,dim_step+1):
      np.savetxt('obsAint_field_step'+str(step)+'.txt', iobs_full[:,:,step])
