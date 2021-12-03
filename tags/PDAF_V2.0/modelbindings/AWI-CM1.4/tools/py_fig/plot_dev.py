import sys
sys.path.append("/home/h/hbkqtang/pyfesom.git/trunk/")
import pyfesom as pf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from matplotlib import cm
from netCDF4 import Dataset
meshpath='/gfs1/work/hbkqtang/input/CORE2_final/'
mesh=pf.load_mesh(meshpath,usepickle=False)
m = Basemap(projection='robin',lon_0=0, resolution='c')
x, y = m(mesh.x2, mesh.y2)

# Read the observation data
fl_obs=Dataset('/gfs1/work/hbkqtang/temporary/obs_SST.nc')

# Read the model forecast from the free run
fl_f=Dataset('/gfs1/work/hbkqtang/output/cpl_work_2014_3/thetao_fesom_20140101.nc')

for i in range(730,1096):
  level_data_obs, elem_no_nan = pf.get_data(fl_obs.variables['sst'][i-730,:],mesh,000)
  level_data_f, elem_no_nan_dummy = pf.get_data(fl_f.variables['thetao'][i,:],mesh,000)

  # Calculate the differnece 
  level_data=level_data_obs-level_data_f

  # Plot the figure
  fig=plt.figure(figsize=(10,7))
  m.drawmapboundary(fill_color='0.9')
  m.drawcoastlines()

  levels = np.arange(-10., 10., 0.05)
  plt.tricontourf(x, y, elem_no_nan[::], level_data, levels = levels, cmap=cm.RdBu_r, extend='both')
  cbar = plt.colorbar(orientation='horizontal', pad=0.03);
  plt.title('SST_obs-fore'+str(i+1))
  plt.tight_layout()
  fig.savefig('SST_obs-fore'+str(i+1))

