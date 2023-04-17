#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt

from math import pi, sin, cos, sqrt, nan

from mpas_tools.io import write_netcdf
#from ap_nearest_gridpoint import ap_nearest_gridpoint
#vector = np.vectorize(np.int_)

#---------------------------------------------
#Get location of data

# E3SM output directory
dir_E3SM_data = '/lustre/scratch4/turquoise/vankova/E3SM/scratch/chicoma-cpu/20230109.GMPAS-JRA1p5-DIB-ISMF.TL319_SOwISC12to60E2r4.chicoma-cpu/run/'

# E3SM restart output file for mesh and coordinates
fname_rst = '20230109.GMPAS-JRA1p5-DIB-ISMF.TL319_SOwISC12to60E2r4.chicoma-cpu.mpaso.rst.0002-01-01_00000.nc'
# E3SM meltrate output file
fname_out_prefix = '20230109.GMPAS-JRA1p5-DIB-ISMF.TL319_SOwISC12to60E2r4.chicoma-cpu.mpaso.hist.am.timeSeriesStatsMonthly.'
fname_out_suffix = '-01.nc'
yr1 = 1 #start year
yrmax = 15 #number of years
jmax = 12 #number of months

dir_fig_save = '/users/vankova/mpas-analysis/iv_analysis/'

# ApRES data directory
#dir_ApRES_data = '~/mpas-analysis/fris_map2/iv_data/ApRES_timeseries/'
dir_ApRES_data = '/usr/projects/e3sm/diagnostics/observations/Ocean/FRIS/ApRES_timeseries/'
# ApRES sites to compare
site_names = ["R02" , "R03" , "R04" , "R05" , "R06" , "R07" , "R08" , "R10" , "R14" , "R15" , "FSW2" , "FSE1" , "FNE3" , "Site5"]
ApRES_prefix = 'FRIS_basal_melt_'
ApRES_loc = f'{dir_ApRES_data}{ApRES_prefix}'

#---------------------------------------------
#---------------------------------------------

#load model coordinates
dsMesh = xarray.open_dataset(f'{dir_E3SM_data}{fname_rst}')
dsMesh = dsMesh[['xCell', 'yCell', 'zCell', 'latCell', 'lonCell']]
dsMesh.load()

x = dsMesh.xCell
y = dsMesh.yCell
z = dsMesh.zCell

def ap_nearest_gridpoint(site_names,dsMesh,ApRES_loc):

    jmax = len(site_names)
    x = dsMesh.xCell
    y = dsMesh.yCell
    z = dsMesh.zCell

    # Find indices of nearest grid point to observations
    ind = np.array(np.zeros((jmax,1)))
    mr_mean = np.array(np.zeros((jmax,1)))
    for j in range(jmax):
        ncApRES_loc = f'{ApRES_loc}{site_names[j]}.nc'
        ncApRES = xarray.open_dataset(ncApRES_loc,decode_times=False)
        ncApRES = ncApRES[['lat', 'lon','mean_melt']]
        ncApRES.load()
        lat = ncApRES.lat
        lon = ncApRES.lon
        mean_melt = ncApRES.mean_melt

        # Compute minimum distance
        erad_e3 = 6.37122e6
        deg2rad = pi / 180

        lon = 360 + lon
        ApRES_z = erad_e3 * sin(lat * deg2rad)
        ApRES_x = erad_e3 * cos(lat * deg2rad) * cos(lon * deg2rad)
        ApRES_y = erad_e3 * cos(lat * deg2rad) * sin(lon * deg2rad)

        #--Distance-----------------------
        dist2point = numpy.absolute(numpy.sqrt((x-ApRES_x)**2 + (y-ApRES_y)**2 + (z-ApRES_z)**2))
        ii = numpy.where(dist2point == dist2point.min())
        ind[j] = np.array(ii)
        mr_mean[j] = np.array(mean_melt)

    return ind, mr_mean

# Find nearest grid point to observations and output index and melt rate value
ind, mr_mean = ap_nearest_gridpoint(site_names, dsMesh, ApRES_loc)

#---------------------------------------------
# Get model timeseries at given indices (nearest lat lon)

kmax = len(ind)
mr = np.array(np.zeros((kmax,jmax*yrmax)))
time = np.array(np.zeros((jmax*yrmax,1)))
ind = ind.astype(int)
ctr = -1;
for y in range(yrmax):
    yr = yr1+y
    for j in range(jmax):
        ctr = ctr + 1;
        #if j + 1 < 10:
        #    fout = f'{dir_E3SM_data}{fname_out_prefix}00{yr}-0{j + 1}{fname_out_suffix}'
        #else:
        #    fout = f'{dir_E3SM_data}{fname_out_prefix}00{yr}-{j + 1}{fname_out_suffix}'
        fout = f'{dir_E3SM_data}{fname_out_prefix}{yr:04}-{j + 1:02}{fname_out_suffix}'
        print(fout)
        out_mr = xarray.open_dataset(fout)
        dsMesh = out_mr[['timeMonthly_avg_landIceFreshwaterFlux']]
        dsMesh.load()
        mr_all = dsMesh.timeMonthly_avg_landIceFreshwaterFlux
        mr_all = np.array(mr_all)
        for k in range(kmax):
            mr[k][ctr] = mr_all[0][ind[k]]/1000*365.25*24*3600
        time[ctr] = (ctr+1)*30-15

# Plot and save modeled vs ApRES melt rate tseries
for k in range(kmax):
    plt.figure()
    ncApRES_loc = f'{ApRES_loc}{site_names[k]}.nc'
    print(ncApRES_loc)
    ncApRES = xarray.open_dataset(ncApRES_loc, decode_times=False)
    ncApRES = ncApRES[['time', 'mean_melt', 'melt_timeseries', 'time_monthly', 'melt_timeseries_monthly']]
    ncApRES.load()
    timeap = ncApRES.time
    t_offset = timeap[0]
    timeap = timeap - t_offset
    timeap_monthly = ncApRES.time_monthly
    timeap_monthly = timeap_monthly - t_offset
    melt = ncApRES.melt_timeseries
    mean_melt = ncApRES.mean_melt
    melt_monthly = ncApRES.melt_timeseries_monthly

    d2y = 365.25
    plt.plot(time/d2y, mr[k][:], label = 'Model')
    #plt.plot(timeap, melt, label = 'ApRES')
    plt.plot(timeap_monthly/d2y, melt_monthly, label='ApRES monthly')
    #plt.plot(timeap, melt*0+mean_melt, label = 'ApRES mean')
    plt.xlabel('Time (years)')
    plt.ylabel('Melt rate (m/a)')
    plt.title(site_names[k])
    plt.legend(loc = 1)
    plt.savefig(f'{dir_fig_save}{site_names[k]}_compare.png', bbox_inches='tight', dpi=300)
    #plt.show()
