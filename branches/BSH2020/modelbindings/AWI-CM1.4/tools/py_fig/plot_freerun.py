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
fl=Dataset('/gfs1/work/hbkqtang/output/cpl_work_2014_3/tos_fesom_20140101.nc')
m = Basemap(projection='robin',lon_0=0, resolution='c')
x, y = m(mesh.x2, mesh.y2)


for i in range(938,1097):
    level_data, elem_no_nan = pf.get_data(fl.variables['tos'][i,:],mesh,000)
    fig=plt.figure(figsize=(10,7))
    m.drawmapboundary(fill_color='0.9')
    m.drawcoastlines()

    levels = np.arange(-5., 35., .1)
    plt.tricontourf(x, y, elem_no_nan[::], level_data, levels = levels, cmap=cm.Spectral_r, extend='both')
    cbar = plt.colorbar(orientation='horizontal', pad=0.03);
    plt.title('SST_free_'+str(i+1))
    plt.tight_layout()
    fig.savefig('SST_free_'+str(i+1))



