

from netCDF4 import Dataset
import numpy as np
import sys
sys.path.append("/home/h/hbkqtang/pyfesom.git/trunk/")
import pyfesom as pf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm

meshpath='/gfs1/work/hbkqtang/input/CORE2_final/'
mesh=pf.load_mesh(meshpath,usepickle=False)

dataset=Dataset('/gfs1/work/hbkqtang/output/cov_2016/cov.nc') 
ens_size=8

nodes_2D=dataset.dimensions['nodes_2D']
nodes_3D=dataset.dimensions['nodes_3D']
n2d=nodes_2D.size
sigma=dataset.variables['sigma'][0:ens_size-1]
sigma=np.float32(sigma)
sigma_trans=np.matrix(sigma)
sigma=np.transpose(sigma_trans)
temp_svd=dataset.variables['temp_svd'][0:ens_size-1,0:n2d]
temp_svd=np.float32(temp_svd)
temp_svd_trans=np.transpose(temp_svd)
cov=np.dot(temp_svd_trans,sigma)
cov=np.dot(cov,sigma_trans)
cov=np.dot(cov,temp_svd)
cov_dia=np.diagonal(cov)
cov_mean=np.mean(cov)
cov_mean_dia=np.mean(cov_dia)
std=np.sqrt(cov_mean_dia)
print(std)

m = Basemap(projection='robin',lon_0=0, resolution='c')
x, y = m(mesh.x2, mesh.y2)

elem2=mesh.elem[mesh.no_cyclic_elem,:]
verbose=True
depth=000


dind=(abs(mesh.zlevs-depth)).argmin()
ind_depth=mesh.n32[:,dind]-1
ind_noempty=np.where(ind_depth>=0)[0]
ind_empty=np.where(ind_depth<0)[0]
level_data=np.zeros(shape=(mesh.n2d))
level_data[ind_noempty]=cov_dia[ind_depth[ind_noempty]]
level_data[ind_empty]=np.nan

d=level_data[elem2].mean(axis=1)
no_nan_triangles = np.invert(np.isnan(d))
elem_no_nan = elem2[no_nan_triangles,:]


fig=plt.figure(figsize=(10,7))
m.drawmapboundary(fill_color='0.9')
m.drawcoastlines()
levels = np.arange(0., .4, .005)
plt.tricontourf(x, y, elem_no_nan[::], level_data, levels = levels, cmap=cm.Spectral_r, extend='both')
cbar = plt.colorbar(orientation='horizontal', pad=0.03);
plt.title('SST_ini_var_hwindow6_4')
plt.tight_layout()
fig.savefig('SST_ini_var_hwindow6_4')

