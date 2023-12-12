from glob import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import cm
from mcradar import getHydroAtmAtt
from cmcrameri import cm as colmap
from attenuationFunctions import * 
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
def getclosestTemp(data,T):
	'''
	calculate range closest to a certain Temperature
	input:
	data: xr.dataset with temperature as variable and range as coordinate
	T: temperature to be selected in Kelvin
	returns: 
	closest range to the input temperature
	'''
	closestT = np.abs(T - data.temperature)
	#allna = closestT.isnull().all(dim='Vel')
	return data.range[closestT.argmin(dim='range')]
def globalAttr(data):
    '''
    adds global Attributes to dataset
    input:
    data: xr.Dataset
    returns:
    data with attributes
    '''
    from datetime import date
    data.attrs['Experiment']= 'TRIPEX-POL-scan, Forschungszentrum Juelich'
    data.attrs['Instrument']= 'joyrad10 (X-Band), joyrad35 (Ka-Band), Wband-scan 94.4 GHz cloud radar (W band)'
    data.attrs['Data']= 'Moments of X, Ka and W-Band, as well as gas attenuation, absolut calibration offset of Ze and attenuation correction of Ze by matching Ze at cloud top (see https://doi.org/10.5194/acp-22-11795-2022 for more explanation)'
    data.attrs['Institution']= 'Data processed within the PROM project FRAGILE, Meteorological Institute, Ludwig-Maximillians University Munich, Germany'
    data.attrs['Latitude']= '50.908547 N'
    data.attrs['Longitude']= '6.413536 E'
    data.attrs['Altitude']= 'Altitude of the JOYCE (www.joyce.cloud) platform: 111m asl'
    data.attrs['process_date'] = str(date.today())
    return data
def generateOffsetDataset(midTime):
	''' 
	generate dataset for offset calculation
	input: 
	midTime: center time of interval in which offset is calculated, pd-daterange object
	returns:
	dataset with offset variables
	''' 
	offsetVar = xr.DataArray(dims=('time','freq'),coords={'time': midTime,'freq':[9.6e9,35.5e9,94e9]})

	dataOffset = xr.Dataset({'offsetMean':offsetVar,
							'offsetStd':offsetVar.copy(),
							'offsetVar':offsetVar.copy(),
							'refRadarStd':offsetVar.copy(),
							'refRadarVar':offsetVar.copy(),
							'validPoints':offsetVar.copy(),
							'totalPoints':offsetVar.copy(),
							'correlation':offsetVar.copy(),
							})
	dataOffset['offsetMean'].attrs = {'description':'mean of offset calculated between X and Ka-Band (select freq=35.5e9) and X and W-Band (select freq=94e9) in 30 min window','units':'dB'}
	dataOffset['offsetStd'].attrs = {'description':'standard deviation of offset calculated between X and Ka-Band (select freq=35.5e9) and X and W-Band (select freq=94e9) in 30 min window','units':'dB'}
	dataOffset['offsetVar'].attrs = {'description':'variance of offset calculated between X and Ka-Band (select freq=35.5e9) and X and W-Band (select freq=94e9) in 30 min window','units':'dB'}
	dataOffset['refRadarStd'].attrs = {'description':'standard deviation of W-Band Ze (reference radar) in 30 min window','units':'dB'}
	dataOffset['refRadarVar'].attrs = {'description':'variance of X-Band Ze (reference radar) in 30 min window','units':'dB'}
	dataOffset['validPoints'].attrs = {'description':'number of valid points in 30 min window','units':''}
	dataOffset['totalPoints'].attrs = {'description':'number of total points in 30 min window','units':''}
	dataOffset['correlation'].attrs = {'description':'correlation between X and Ka-Band Ze (select freq=35.5e9) and between X and W-band (select freq=94e9) 30 min window','units':''}
	
	return dataOffset						
def getInterval(data,interval=[-30,-10]):
	'''
	get correct interval for Ze matching, 1km above melting layer and the interval defined in interval
	input:
	data: xr.dataset containing Ze with freq as coordinate and temperature as variable
	interval: interval in which Ze are used 
	returns:
	data which is masked everwhere outside interval and below 1km above ML
	'''
	#- get range of 1km above ML:
	MLrange = getclosestTemp(data,273.15)
	dataCold = data.where(data.range > MLrange+1000)
	
	#- select only cases where Ze within predifend interval: -30 to -10dB for Ka-W (Ka-Band reflectivity), -15 to 0 for X-Ka (Ka-Band reflectivity)
	dataInt = dataCold.where(dataCold.Ze.sel(freq=9.6e9)>interval[0])
	dataInt = dataInt.where(dataInt.Ze.sel(freq=9.6e9)<interval[1])

	return dataInt

def getTime(date):
	'''
	get time interval arrays in which to make the matching:
	input:
	date: pd.datetime 
	returns:
	beginTime: starttime of interval
	midTime: middle point of interval
	endTime: endtime of interval
	'''
	#TODO: this should be more flexible dealing with interval size!!
	beginTime = pd.date_range(date,date+pd.offsets.Day(1)-pd.offsets.Minute(30)-pd.offsets.Second(4),freq='1T')
	midTime = pd.date_range(date + pd.offsets.Minute(15),date+pd.offsets.Day(1)-pd.offsets.Minute(15),freq='1T')
	endTime = pd.date_range(date + pd.offsets.Minute(30),date+pd.offsets.Day(1)-pd.offsets.Second(4),freq='1T')
	#15min
	#beginTime = pd.date_range(date,date+pd.offsets.Day(1)-pd.offsets.Minute(15)-pd.offsets.Second(4),freq='5T')
	#midTime = pd.date_range(date + pd.offsets.Minute(7)+pd.offsets.Second(30),date+pd.offsets.Day(1)-pd.offsets.Minute(7)-pd.offsets.Second(34),freq='5T')
	#endTime = pd.date_range(date + pd.offsets.Minute(15),date+pd.offsets.Day(1)-pd.offsets.Second(4),freq='5T')
	#45 min
	#beginTime = pd.date_range(date,date+pd.offsets.Day(1)-pd.offsets.Minute(45)-pd.offsets.Second(4),freq='5T')
	#midTime = pd.date_range(date + pd.offsets.Minute(22)+pd.offsets.Second(30),date+pd.offsets.Day(1)-pd.offsets.Minute(22)-pd.offsets.Second(34),freq='5T')
	#endTime = pd.date_range(date + pd.offsets.Minute(45),date+pd.offsets.Day(1)-pd.offsets.Second(4),freq='5T')
	#60min
	#beginTime = pd.date_range(date,date+pd.offsets.Day(1)-pd.offsets.Minute(60)-pd.offsets.Second(4),freq='5T')
	#midTime = pd.date_range(date + pd.offsets.Minute(30),date+pd.offsets.Day(1)-pd.offsets.Minute(30),freq='5T')
	#endTime = pd.date_range(date + pd.offsets.Minute(60),date+pd.offsets.Day(1)-pd.offsets.Second(4),freq='5T')
	#2h
	#beginTime = pd.date_range(date,date+pd.offsets.Day(1)-pd.offsets.Minute(120)-pd.offsets.Second(4),freq='1T')
	#midTime = pd.date_range(date + pd.offsets.Minute(60),date+pd.offsets.Day(1)-pd.offsets.Minute(30),freq='1T')
	#endTime = pd.date_range(date + pd.offsets.Minute(120),date+pd.offsets.Day(1)-pd.offsets.Second(4),freq='1T')

	return beginTime,midTime,endTime
def getOffset(dataXKa,dataXW,date):
	'''
	calculates offset between Ze(Ka-Band) and Ze (x or W-Band)
	input:
	dataXKa: data in correct interval for XKa matching(from getInterval)
	dataXW: data in correct interval for XW matching (from getInterval)
	returns:
	dataOffset: dataset containing offsets
	'''
	#- get time interval
	beginTime,midTime,endTime=getTime(date)
	
	#- generate offset dataset:	
	dataOffset = generateOffsetDataset(midTime)
	
	for i,tstart,tmid,tend in zip(range(len(beginTime)),beginTime,midTime,endTime):
		ZeXWsel = dataXW.Ze.sel(time=slice(tstart,tend))
		ZeXKasel = dataXKa.Ze.sel(time=slice(tstart,tend))
		offsetW = ZeXWsel.sel(freq=9.6e9) - ZeXWsel.sel(freq=94e9)
		offsetKa = ZeXKasel.sel(freq=9.6e9) - ZeXKasel.sel(freq=35.5e9)
		dataOffset['offsetMean'].loc[tmid,94e9] = offsetW.mean()
		dataOffset['offsetStd'].loc[tmid,94e9] = offsetW.std()
		dataOffset['offsetVar'].loc[tmid,94e9] = offsetW.var()
		dataOffset['validPoints'].loc[tmid,94e9] = offsetW.count()
		dataOffset['totalPoints'].loc[tmid,94e9] = len(offsetW.values.flatten())
		dataOffset['correlation'].loc[tmid,94e9] = xr.corr(ZeXWsel.sel(freq=9.6e9),ZeXWsel.sel(freq=94e9))
		dataOffset['refRadarStd'].loc[tmid,94e9] = ZeXWsel.sel(freq=9.6e9).std()
		
		dataOffset['offsetMean'].loc[tmid,35.5e9] = offsetKa.mean()
		dataOffset['offsetStd'].loc[tmid,35.5e9] = offsetKa.std()
		dataOffset['offsetVar'].loc[tmid,35.5e9] = offsetKa.var()
		dataOffset['validPoints'].loc[tmid,35.5e9] = offsetKa.count()
		dataOffset['totalPoints'].loc[tmid,35.5e9] = len(offsetKa.values.flatten())
		dataOffset['correlation'].loc[tmid,35.5e9] = xr.corr(ZeXKasel.sel(freq=9.6e9),ZeXKasel.sel(freq=35.5e9))
		dataOffset['refRadarStd'].loc[tmid,35.5e9] = ZeXKasel.sel(freq=9.6e9).std()
		
		dataOffset['offsetMean'].loc[tmid,9.6e9] = 0
		dataOffset['offsetStd'].loc[tmid,9.6e9] = 0
		dataOffset['offsetVar'].loc[tmid,9.6e9] = 0
		dataOffset['validPoints'].loc[tmid,9.6e9] = 0
		dataOffset['totalPoints'].loc[tmid,9.6e9] = 0
		dataOffset['correlation'].loc[tmid,9.6e9] = 0
		dataOffset['refRadarStd'].loc[tmid,9.6e9] = 0
		
		print(tend)
	
	return dataOffset

#def getVarianceFlag(dataXKa,dataKaW):
#	
#	offsetVar = dataXKa.Ze.copy()*np.nan
#	print(offsetVar)
#	offsetVar.loc[9.6e9,:,:] = (dataXKa.Ze.sel(freq=9.6e9)-dataXKa.Ze.sel(freq=35.5e9)).rolling(time=int(np.floor(30*60/4)), min_periods=1, center=True).var()
#	offsetVar.loc[94e9,:,:] = (dataKaW.Ze.sel(freq=35.5e9)-dataKaW.Ze.sel(freq=94e9)).rolling(time=int(np.floor(30*60/4)), min_periods=1, center=True).var()
#	#data['varFlag'] = data.varFlag.fillna(3)
#	varFlag = (offsetVar>2).any(dim='range')
#	varFlag.attrs = {'description':'flag where the variance of difference in any rangegate is larger than 2. True if larger any are larger than 2.'}
#	offsetVar.attrs = {'description':'variance of offset calculated in 30 min rolling window for X and Ka-Band (select freq=9.6e9) and W and Ka-Band (select freq=94e9)','units':'dB2'}
#	
#	return offsetVar,varFlag

def getVarianceFlag(dataXKa,dataKaW,fillna=0):
	'''
	get variance in time of offset (DWR) to potentially masked areas where variance is too large (variance flag in Joses retrieval)
	input: 
	data: xr.dataset with Ze as variable
	returns:
	data: dataset with added variance and varFlag
	'''
	offsetVar = dataKaW.Ze.copy()*np.nan
	print(offsetVar)
	offsetVar.loc[35.5e9,:,:] = (dataXKa.Ze.sel(freq=9.6e9)-dataXKa.Ze.sel(freq=35.5e9)).rolling(time=int(np.floor(30*60/4)), min_periods=1, center=True).var()
	offsetVar.loc[94e9,:,:] = (dataKaW.Ze.sel(freq=9.6e9)-dataKaW.Ze.sel(freq=94e9)).rolling(time=int(np.floor(30*60/4)), min_periods=1, center=True).var()
	offsetVar = offsetVar.fillna(fillna)
	#data['varFlag'] = data.varFlag.fillna(3)
	varFlag = offsetVar.where(offsetVar > 2,0)# make everything that is smaller than 2 equal to 0 #(offsetVar>2).any(dim='range')
	varFlag = varFlag.where(offsetVar <= 2,1)# make everything larger than 2 equal to 1
	#varFlag = varFlag.astype(bool)
	
	#fig,ax = plt.subplots(nrows=2)
	#cmap = ListedColormap(["blue", "green","red"])
	#p1=ax[0].pcolormesh(offsetVar.time,offsetVar.range,offsetVar.sel(freq=94e9).T,cmap=cmap,vmin=0,vmax=3)
	#plt.colorbar(p1,ax=ax[0])
	#cmap = ListedColormap(["blue", "green"])
	#p1=ax[1].pcolormesh(varFlag.time,varFlag.range,varFlag.sel(freq=94e9).T,cmap=cmap,vmin=-0.5,vmax=1.5)
	#plt.colorbar(p1,ax=ax[1])
	#plt.show()
	
	varFlag.attrs = {'description':'flag where the variance of difference is larger than 2. True if larger than 2.'}
	offsetVar.attrs = {'description':'variance of offset calculated in 30 min rolling window for X and Ka-Band (select freq=35.5e9) and X- and W-Band (select freq=94e9)','units':'dB2'}
	#data['varRolling'] = offsetVar
	#data['varFlag'] = varFlag
	return varFlag#offsetVar,varFlag

def getCategory(variable,flag,bitvalue,multiple=False):
# get categories of CN. with bitvalue you can specify which category(ies) you want: 
# 0: cloud droplets, 1:falling hydrometeors, 2:colder then 0°C, 3: melting ice particles, 4: aerosols, 5: insects
# by giving multiple bitvalues like [0,1,2] and adding multiple=True, you can get areas where multiple conditions apply. 
# e.g. with [0,1,2] you can get cloud regions where we have cloud droplets along with freezing hydrometeors (so ice)
# you need to specify a variable where the mask is applied to. most simply would be the radar reflectivity, because this usually has the shape of the clouds
  #((flag.values >> bitvalue) & 1).astype(bool)
  #print(flag)
  #print(bitvalue)
  #variable = variable.where(((flag.values >> 3) & 1).astype(bool))
  #except:
  #  print(variable,' has time problem')  
  #return variable
	mask = int('000000',2)
	if multiple == True:
		for bit in bitvalue:
			mask |= (1<<bit)
	else:
		bit = bitvalue
		mask |= (1<<bit)
	#try:  
	variable = variable.where((flag.values & mask) == mask)
	#except:
	#  print(variable,' has time problem')  
	return variable  

def maskData(variable,flag,bit):
	return variable.where(~((flag.values >> bit) & 1).astype(bool))
	
def getFlags(data):
	'''
	this gets the quality flags as in Jose: 
	bit 6: LWP of CN larger than 200gm-2;
	bit 7: rain detected by cloudnet
	bit 13: variance in time of DWR>2dB
	bit 14: correlation of points is smaller than 0.7
	bit 15: number of valid measurements is smaller than 300
	INPUT:
	data: xr.dataset after attenuation correction was done (after function getOffset)
	RETURNS:
	data: xr.dataset containing the bit flag 
	'''
	#- get variance flag, only in region where we calculated the offset!
	dataKaW = getInterval(data,interval=[-30,-10])
	dataXKa = getInterval(data,interval=[-15,0])
	varFlag = getVarianceFlag(dataXKa,dataKaW)
	#varFlag2 = getVarianceFlag(dataXKa,dataKaW,fillna=3)
	#- get correlation flag
	corrFlag = data.correlation.where(data.correlation <= 0.7,1) #everything that is larger than 0.7 will be 1 (so True to be flagged)
	corrFlag = corrFlag.where(data.correlation > 0.7,0) # everything that is smaller than 0.7 will be 0 (so false to be flagged)
	corrFlag = -1*(corrFlag-1) # need to switch 0 and 1 because we want to have everything above 0.7 to not be flagged
	corrFlag = corrFlag.expand_dims({'range':data.range})#.transpose('freq','time','range') # this is necessary because variance flag has range as dim!
	#- get number of valid points flag
	numFlag = data.validPoints.where(data.validPoints <= 300,1) #everything that is larger than 300 will be 1 (so True to be flagged)
	numFlag = numFlag.where(data.validPoints > 300,0)# everything that is smaller than 300 will be 0 (so false to be flagged)
	numFlag = -1*(numFlag-1) # need to switch 0 and 1 because we want to have everything above 0.7 to not be flagged
	numFlag = numFlag.expand_dims({'range':data.range})#.transpose('freq','time','range')
	#- now lets get rainflag:
	#- get cloudnet categories and select only liquid precipitation
	drizzle_rain = getCategory(data.category_bits/data.category_bits,data.category_bits.astype('uint16'),1)
	drizzle_rain_drop = getCategory(data.category_bits/data.category_bits,data.category_bits.astype('uint16'),[0,1],multiple=True)
	ice = getCategory(data.category_bits/data.category_bits,data.category_bits.astype('uint16'),[1,2],multiple=True)
	drizzle_rain = drizzle_rain.where(ice!=1)
	drizzle_rain_drop = drizzle_rain_drop.where(ice!=1)
	rain = drizzle_rain.fillna(0) + drizzle_rain_drop.fillna(0)
	rain = rain.where(rain<1,1)
	rainFlag = rain.where(rain==1,np.nan).any(dim='range') # have flag when there is rain detected in any range gate
	rainFlag = rainFlag.expand_dims({'range':data.range,'freq':data.freq})#.transpose('freq','time','range')
	#- get lwp flag from CN, flag st to true when lwp > 200
	lwpFlag = data.lwp.where(data.lwp <= 200,1)# make everything that is larger than 200 equal to 1 #(offsetVar>2).any(dim='range')
	lwpFlag = lwpFlag.where(data.lwp > 200,0)# make everything smaller than 200 equal to 0
	lwpFlag = lwpFlag.expand_dims({'range':data.range,'freq':data.freq})#.transpose('freq','time','range')
	#- now make it into bits:
	lwpFlagBit = np.array(lwpFlag.transpose('freq','time','range').values.astype('uint16'))<<6
	rainFlagBit = np.array(rainFlag.transpose('freq','time','range').values.astype('uint16'))<<7
	varFlagBit = np.array(varFlag.transpose('freq','time','range').values.astype('uint16'))<<13
	corrFlagBit = np.array(corrFlag.transpose('freq','time','range').values.astype('uint16'))<<14
	numFlagBit = np.array(numFlag.transpose('freq','time','range').values.astype('uint16'))<<15
	# add all flags together and put it back into an xarray dataarray
	finalFlag = rainFlagBit + lwpFlagBit + varFlagBit + corrFlagBit + numFlagBit 
	finalFlag = xr.DataArray(finalFlag,dims=('freq','time','range'),coords={'freq':data.freq,'time':data.time,'range':data.range})
	finalFlag.attrs = {'description':'quality flag indicating reliability of the offset correction. ',
						'comments':'For W-Band offset flags select freq=94e9, for Ka-Band offset flags select freq=35.5e9. Bits 0 to 5: empty; Bit 6: 0 if the liquid water path is less than 200 g m-2, 1 if the liquid water path is greater than 200 g m-2; Bit 7: 0 no rain, 1 rain; Bits 8 to 12: empty; Bit 13: 0 if the variance of the DWR_X_Ka within a 30 minutes time window is less than 2 dB**2, 1 if the DWR variance is greater than 2 dB**2; Bit 14: 0 if correlation of X and Ka/W Band reflectivities is larger than 0.7, 1 if correlation is less than 0.7; Bit 15: 0 if the number of points used to calculate the offset is greater than 300, 1 if the number of points is less than 300. If Bit 14 or higher is set, we recommend not to use the calculated offsets but e.g. rather interpolate between time periods with high-quality offset estimates.'}
	#data['offset_flag'] = finalFlag
	'''
	fig,ax=plt.subplots(nrows=5,figsize=(20,25),sharex=True)
	var = maskData(data.Ze.sel(freq=9.6e9),data.offset_flag.sel(freq=9.6e9),6)
	p1=ax[0].pcolormesh(data.time,data.range,var.T,vmin=-30,vmax=20,cmap=colmap.batlow)
	ax[0].plot(data.time,lwpFlag.where(lwpFlag==1).sel(range=36,freq=9.6e9,method='nearest')*10000,c='C0',lw=3)
	cb=plt.colorbar(p1,ax=ax[0],pad=0.01,aspect=30)
	cb.set_label('Ze X-Band [dBz]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	ax[0].set_title('lwp-flag',fontsize=24)
	
	var = maskData(data.Ze.sel(freq=9.6e9),data.offset_flag.sel(freq=9.6e9),7)
	p1=ax[1].pcolormesh(data.time,data.range,var.T,vmin=-30,vmax=20,cmap=colmap.batlow)
	ax[1].plot(data.time,rainFlag.where(rainFlag==1).sel(range=36,freq=9.6e9,method='nearest')*10000,c='C0',lw=3)
	cb=plt.colorbar(p1,ax=ax[1],pad=0.01,aspect=30)
	cb.set_label('Ze X-Band [dBz]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	ax[1].set_title('rain-flag',fontsize=24)
	
	var = maskData(data.Ze.sel(freq=9.6e9),data.offset_flag.sel(freq=9.6e9),13)
	p1=ax[2].pcolormesh(data.time,data.range,var.T,vmin=-30,vmax=20,cmap=colmap.batlow)
	#ax[2].plot(data.time,varFlag.where(lwpFlag==1).sel(range=36,freq=9.6e9,method='nearest')*10000,c='C0')
	cb=plt.colorbar(p1,ax=ax[2],pad=0.01,aspect=30)
	cb.set_label('Ze X-Band [dBz]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	ax[2].set_title('Var-flag',fontsize=24)
	
	var = maskData(data.Ze.sel(freq=9.6e9),data.offset_flag.sel(freq=9.6e9),14)
	p1=ax[3].pcolormesh(data.time,data.range,var.T,vmin=-30,vmax=20,cmap=colmap.batlow)
	ax[3].plot(data.time,corrFlag.where(corrFlag==1).sel(range=36,freq=94e9,method='nearest')*10000,c='C0',lw=3)
	cb=plt.colorbar(p1,ax=ax[3],pad=0.01,aspect=30)
	cb.set_label('Ze X-Band [dBz]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	ax[3].set_title('Corr-flag',fontsize=24)
	
	var = maskData(data.Ze.sel(freq=9.6e9),data.offset_flag.sel(freq=9.6e9),15)
	p1=ax[4].pcolormesh(data.time,data.range,var.T,vmin=-30,vmax=20,cmap=colmap.batlow)
	ax[4].plot(data.time,numFlag.where(numFlag==1).sel(range=36,freq=94e9,method='nearest')*10000,c='C0',lw=3)
	cb=plt.colorbar(p1,ax=ax[4],pad=0.01,aspect=30)
	cb.set_label('Ze X-Band [dBz]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	ax[4].set_title('Number-flag',fontsize=24)
	for i,a in enumerate(ax):
		a.tick_params(labelsize=18)
		a.set_ylabel('range [m]',fontsize=20)
	plt.tight_layout()
	plt.savefig('plots/20220206_test_bitflags_X-band.png')
	plt.close()
	print(finalFlag)
	quit()
	'''
	return finalFlag
			
def getJoyrad(pathData,date,rangeRef,timeRef,band,offset=0):
	'''
	get Ka or X Band data with correct shape
	input:
	pathData: path where data is
	date: date to process (pd. datetime)
	rangeRef: reference of range to reindex
	timeRef: reference of time (pd daterange)
	band: X or Ka
	offset: calibration offset calculated with raincoat
	returns:
	data: dataset containing radar moments and Ze with calibration offset
	'''
	if band == 'Ka':
		joyrad='35'
		freq = 35.5e9
	elif band == 'X':
		joyrad='10'
		freq = 9.6e9
	data = xr.open_dataset('{path}campaign_aux_data/tripex-pol/joyrad{joyrad}/resampled/{date}_mon_joyrad{joyrad}.nc'.format(path=pathData,date=date.strftime('%Y%m%d'),joyrad=joyrad))
	print(data)
	data = data.reindex({'range':rangeRef},method='nearest',tolerance=15)
	data = data.reindex({'time':timeRef},method='nearest',tolerance='4S').expand_dims({'freq':np.array([freq])}).rename({'Zg':'Ze','VELg':'MDV','RMSg':'WIDTH','SKWg':'SK'})
	data['calibration_offset'] = xr.DataArray(data=offset,dims=['freq'],coords={'freq':np.array([freq])})
	data['calibration_offset'].attrs = {'description':'absolute calibration offset obtained with raincoat by comparing radar Ze with forward simulated Ze from Parsivel DSD measurements','units':'dB','comment':'positive means radar is smaller than parsivel'}
	
	return data
	
def getWband(pathData,date,offset=0):
	'''
	get Ka or X Band data with correct shape
	input:
	pathData: path where data is
	date: date to process (pd. datetime)
	rangeRef: reference of range to reindex
	timeRef: reference of time (pd daterange)
	band: X or Ka
	offset: calibration offset calculated with raincoat
	returns:
	data: dataset containing radar moments and Ze with calibration offset
	'''
	
	data = xr.open_dataset('{path}{date}_ZEN_moments_wband_scan.nc'.format(path=pathData,date=date.strftime('%Y%m%d')))
	print(data)
	data = data[['Ze','MDV','WIDTH','SK']].expand_dims({'freq':np.array([94e9])}) 
	data['calibration_offset'] = xr.DataArray(data=offset,dims=['freq'],coords={'freq':np.array([94e9])})
	data['calibration_offset'].attrs = {'description':'absolute calibration offset obtained with raincoat by comparing radar Ze with forward simulated Ze from Parsivel DSD measurements','units':'dB','comment':'positive means radar is smaller than parsivel'}
	
	return data
