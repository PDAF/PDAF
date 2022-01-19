

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

dataset=Dataset('/home/h/hbkqtang/fesom_echam6_oasis3-mct/run_scripts/hlrn3/work/cpl_work_8sst_testcase/fesom.2016.oce.daexp1.nc') 
temp_std=dataset.variables['temp_ini'][0:1,0:126859]
print(temp_std.min())
print(temp_std.max())
temp_std_mean=np.mean(temp_std)
temp_std_std=np.sqrt(temp_std_mean)
print(temp_std_std)

m = Basemap(projection='robin',lon_0=0, resolution='c')
x, y = m(mesh.x2, mesh.y2)

level_data, elem_no_nan = pf.get_data(dataset.variables['temp_ini'][0,:],mesh,000)
fig=plt.figure(figsize=(10,7))
m.drawmapboundary(fill_color='0.9')
m.drawcoastlines()

levels = np.arange(0., 2., .05)
plt.tricontourf(x, y, elem_no_nan[::], level_data, levels = levels, cmap=cm.Spectral_r, extend='both')
cbar = plt.colorbar(orientation='horizontal', pad=0.03);
plt.title('SST_var_model_2')
plt.tight_layout()
fig.savefig('SST_var_model_2')

