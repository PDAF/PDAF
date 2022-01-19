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
elem2=mesh.elem[mesh.no_cyclic_elem,:]
verbose=True
depth=000
dind=(abs(mesh.zlevs-depth)).argmin()
ind_depth=mesh.n32[:,dind]-1
ind_noempty=np.where(ind_depth>=0)[0]
ind_empty=np.where(ind_depth<0)[0]




# Read the observation data
fl_obs=Dataset('/gfs1/work/hbkqtang/temporary/obs_SST.nc')

# Read the model forecast from the free run
fl_f=Dataset('/gfs1/work/hbkqtang/output/cpl_work_2016_1year/thetao_fesom_20160101.nc')

nodes_2D=fl_obs.dimensions['nodes_2D']
n2d=nodes_2D.size
level_obs=np.zeros(shape=(mesh.n2d))
level_fore=np.zeros(shape=(mesh.n2d))

for i in range(0,366):
  if i % 5 == 0:
    obs_data=fl_obs.variables['sst'][i:i+5,:]
    obs_mean=np.nanmean(obs_data,axis=0)
    level_obs[ind_noempty]=obs_mean[ind_depth[ind_noempty]]
    level_obs[ind_empty]=np.nan
    d=level_obs[elem2].mean(axis=1)
    no_nan_triangles = np.invert(np.isnan(d))
    elem_no_nan = elem2[no_nan_triangles,:]
    fore_data=fl_f.variables['thetao'][i:i+5,0:n2d]
    fore_mean=np.nanmean(fore_data,axis=0)
    level_fore[ind_noempty]=fore_mean[ind_depth[ind_noempty]]
    level_fore[ind_empty]=np.nan
    # Calculate the differnece 
    level_data=obs_mean-fore_mean
    # Plot the figure
    fig=plt.figure(figsize=(10,7))
    m.drawmapboundary(fill_color='0.9')
    m.drawcoastlines()
    levels = np.arange(-2., 2., .01)
    plt.tricontourf(x, y, elem_no_nan[::], level_data, levels = levels, cmap=cm.RdBu_r, extend='both')
    cbar = plt.colorbar(orientation='horizontal', pad=0.03);
    plt.title('SST_5day_obs-fore'+str(i+1))
    plt.tight_layout()
    fig.savefig('SST_5day_obs-fore'+str(i+1))
  
