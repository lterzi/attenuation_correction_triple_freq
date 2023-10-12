#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
import xarray as xr
import pandas as pd
def maskData(variable, flag, mask):
    
    variable = variable.where((flag.values & mask) != mask)
    
    return variable
def plot_DWR_offset(data, plotOutPath, strDate, plotID):
    DWR_X_Ka = data['X_DBZ_H'] - data['Ka_DBZ_H']
    DWR_Ka_W = data['Ka_DBZ_H'] - data['W_DBZ_H']

    qFlagW = data['quality_flag_offset_w'].copy()
    qFlagX = data['quality_flag_offset_x'].copy()
    offsetW = data['offset_w'].copy()
    offsetX = data['offset_x'].copy()

    selFlagX = qFlagX.sel(range=4000,method='Nearest')
    selFlagW = qFlagW.sel(range=4000,method='Nearest')

    selOffsetW = offsetW.sel(range=4000,method='Nearest')
    selOffsetX = offsetX.sel(range=4000,method='Nearest')

    fig, axes = plt.subplots(nrows=3, figsize=(18,12),sharex=True)
    radData = {'DWR_X_Ka':{'data':DWR_X_Ka, 'axis':axes[0], 'lim':(-5,20)},
               'DWR_Ka_W':{'data':DWR_Ka_W, 'axis':axes[1], 'lim':(-5,20)}}
    for rad in radData.keys():

        plot = radData[rad]['data'].T.plot(ax=radData[rad]['axis'],
                                        vmax=radData[rad]['lim'][1],
                                        vmin=radData[rad]['lim'][0],
                                        cmap='jet',label='_nolegend_')

        radData[rad]['axis'].set_title(rad +' '+ strDate)
        plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
        radData[rad]['axis'].grid()
        radData[rad]['axis'].set_xlabel('')
        radData[rad]['axis'].set_ylabel('height in km')
    
    
    #colors = [(0,0,0),(0,0,0)]
    #cmap_name = 'fake_white'
    #newcmp = LinearSegmentedColormap.from_list(cmap_name, colors, N=2) 
# fake up the array of the scalar mappable. Urgh…
    #newcmp._A = []
    cmap = ListedColormap(["darkorange", "gold", "lawngreen", "lightseagreen"])

    selOffsetW.plot(color='Orange',label='Offset W',lw=2)
    selOffsetX.plot(color='C0',label='Offset X',lw=2)
    selOffsetW.where(selFlagW<16000).plot(color='Orange',lw=6)
    selOffsetX.where(selFlagX<16000).plot(color='C0',lw=6)
    axes[2].set_ylabel('height in km')
    axes[2].set_title('')
    plt.colorbar(cmap)
    plt.legend()
    plotFileName = ('_').join([strDate,plotID])
    filePathName = ('/').join([plotOutPath,plotFileName])
    plt.show()
    #plt.savefig(filePathName+'.png',dpi=200, bbox_inches='tight')
############################################################################################################################

#dataPath = '/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_2/'
fileID = '_tripex_pol_3fr_L2_mom.nc'
outputPath = '/work/lvonterz/tripexProcessing/output/'
dataPath = '/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_2/' 
dateStart = pd.to_datetime('20181101'); dateEnd = pd.to_datetime('20181231')
dateList = pd.date_range(dateStart, dateEnd,freq='D')
for date in dateList:
    dateID = date.strftime('%Y%m%d')
    data = xr.open_dataset(dataPath+dateID+fileID)
    plot_DWR_offset(data,outputPath, dateID, 'LVL2_DWR_offset')

'''
    qFlagW = data['quality_flag_offset_w'].copy()
    qFlagX = data['quality_flag_offset_x'].copy()
    offsetW = data['offset_w'].copy()
    offsetX = data['offset_x'].copy()

    selFlagX = qFlagX.sel(range=4000,method='Nearest')
    selFlagW = qFlagW.sel(range=4000,method='Nearest')

    selOffsetW = offsetW.sel(range=4000,method='Nearest')
    selOffsetX = offsetX.sel(range=4000,method='Nearest')
    print(offsetW)
    selOffsetW.plot()
    plt.show()
    quit()
   
    selOffsetW.plot(color='Orange',label='Offset W',lw=2)
    selOffsetX.plot(color='C0',label='Offset X',lw=2)
    selOffsetW.where(selFlagW<16000).plot(color='Orange',lw=6)
    selOffsetX.where(selFlagX<16000).plot(color='C0',lw=6)
'''
    







