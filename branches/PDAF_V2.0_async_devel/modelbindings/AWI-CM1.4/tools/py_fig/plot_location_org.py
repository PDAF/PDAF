import sys
sys.path.append("/home/h/hbkqtang/pyfesom.git/trunk/")
import pyfesom as pf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from matplotlib import cm
from netCDF4 import Dataset
m = Basemap(projection='robin',lon_0=0, resolution='c')
en4=Dataset("/gfs1/work/hbkqtang/temporary/EN.4.2.1.f.profiles.g10.201601.nc",format="NETCDF3_64BIT_OFFSET")
# find out the number of profiles for one day
n_prof=en4.dimensions['N_PROF']
n_prof=n_prof.size
lon=en4.variables['LONGITUDE'][:]
lat=en4.variables['LATITUDE'][:]
depth=en4.variables['DEPH_CORRECTED'][:,:]
temp=en4.variables['PSAL_CORRECTED'][:,:]
qc=en4.variables['PROFILE_PSAL_QC']
judate=en4.variables['JULD'][:]
refjud=en4.variables['REFERENCE_DATE_TIME']
Y=refjud[0:4]
Y = [int(x) for x in Y ]
Y=Y[0]*1000+Y[1]*100+Y[2]*10+Y[3]
M=refjud[4:6]
M = [int(x) for x in M ]
M=M[0]*10+M[1]
D=refjud[6:8]
D = [int(x) for x in D ]
D=D[0]*10+D[1]
HOUR=refjud[8:10]
HOUR = [int(x) for x in HOUR ]
HOUR=HOUR[0]*10+HOUR[1]
MINUTE=refjud[10:12]
MINUTE = [int(x) for x in MINUTE ]
MINUTE=MINUTE[0]*10+MINUTE[1]
SECOND=refjud[12:14]
SECOND = [int(x) for x in SECOND ]
SECOND=SECOND[0]*10+SECOND[1]
refjudate=2433283
for date in range(3,4):
     nprof_day=0
     index_obs=np.empty(1500)
     for i in range(n_prof):
        JD=judate[i]+refjudate
        L=int(JD)+68569
        N=int(4*L/146097)
        L=int(L-(146097*N+3)/4)
        Y=int(4000*(L+1)/1461001)
        L=int(L-1461*Y/4+31)
        M=int(80*L/2447)
        D=int(L-2447*M/80)
        L=int(M/11)
        M=int(M+2-12*L)
        Y=int(100*(N-49)+Y+L)
        YEAR=Y
        MONTH=M
        DAY=D
        if (DAY==date and int(qc[i]) ==1):
        #if DAY==date:
               nprof_day=nprof_day+1
               index_obs[nprof_day-1]=i
     lon_day=np.empty(nprof_day)
     lat_day=np.empty(nprof_day)
     depth_day=np.empty([nprof_day,400])
     temp_day=np.empty([nprof_day,400])
     for j in range(nprof_day):
        lon_day[j]=lon[int(index_obs[j])]
        lat_day[j]=lat[int(index_obs[j])]
        depth_day[j,:]=depth[int(index_obs[j]),:]
        temp_day[j,:]=temp[int(index_obs[j]),:]
     temp_100=99999*np.ones(nprof_day)
     #temp_100=99999.0
     for k in range(nprof_day):
        for mm in range (400):
           if (np.abs((depth_day[k,mm]-100.0)))<0.5:
                temp_100[k]=temp_day[k,mm]
     cnt=0
     for nn in range(nprof_day):
        if temp_100[nn]<50.0:
           cnt=cnt+1
     temp_real=np.empty(cnt)
     lon_real=np.empty(cnt)
     lat_real=np.empty(cnt)
     cnt=0
     for nn in range(nprof_day):
        if temp_100[nn]<50.0:
          cnt=cnt+1
          temp_real[cnt-1]=temp_100[nn]
          lon_real[cnt-1]=lon_day[nn]
          lat_real[cnt-1]=lat_day[nn]

x,y=m(lon_real,lat_real)
m.scatter(x,y,c=temp_real,marker='o',edgecolors='none',s=20,cmap=cm.Spectral_r)
plt.clim(31.0, 38.0)
plt.colorbar(orientation='vertical', pad=0.03)
plt.tight_layout()
plt.show()

