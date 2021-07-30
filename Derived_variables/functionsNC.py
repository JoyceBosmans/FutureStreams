import numpy as np
import netCDF4 as nc
import os
import datetime

# Created by Niko Wanders (n.wanders@uu.nl)

MV = 1e20
smallNumber = 1E-39

def getNCData(ncFile, varName, date=None):    
  f = nc.Dataset(ncFile)
  if date != None:
    nctime = f.variables['time']  # A netCDF time variable object.
    idx = nc.date2index(date, nctime, calendar = nctime.calendar, \
                                                  select='exact')
    data = f.variables[varName][idx,:,:]
  else:
    data = f.variables[varName][:]
  f.close()
  return data

def getNCIndexData(ncFile, varName = 'discharge', idxStart=None, idxEnd=None):
  f = nc.Dataset(ncFile)
  if idxStart != None:
    data = f.variables[varName][idxStart:idxEnd,:,:]
  else:
    data = f.variables[varName][:]
  f.close()
  return data

def readTimeValues(ncFile):
  # Get netCDF file and variable name:
  f = nc.Dataset(ncFile)
  nctime = f.variables['time'][:]
  nctimeUnit = f.variables['time'].units
  nctimeCalendar = f.variables['time'].calendar
  timeVar = nc.num2date(nctime,units = nctimeUnit, calendar = nctimeCalendar)
  return timeVar
  
def createNetCDF(ncFileName, varName, varUnits, latitudes, longitudes,\
                                      longName = None, loop=False, title="Bioclimoutput"):
    
    rootgrp= nc.Dataset(ncFileName,'w', format='NETCDF4')
    
    #-create dimensions - time is unlimited, others are fixed
    rootgrp.createDimension('time',None)
    rootgrp.createDimension('lat',len(latitudes))
    rootgrp.createDimension('lon',len(longitudes))
    
    date_time= rootgrp.createVariable('time','f4',('time',))
    date_time.standard_name= 'time'
    date_time.long_name= 'Days since 1901-01-01'
    
    date_time.units= 'Days since 1901-01-01' 
    date_time.calendar= 'standard'
        
    lat= rootgrp.createVariable('lat','f4',('lat',))
    lat.long_name= 'latitude'
    lat.units= 'degrees_north'
    lat.standard_name = 'latitude'
    
    lon= rootgrp.createVariable('lon','f4',('lon',))
    lon.standard_name= 'longitude'
    lon.long_name= 'longitude'
    lon.units= 'degrees_east'
    
    lat[:]= latitudes
    lon[:]= longitudes
    
    if loop:
        for i in range(len(varName)):
            shortVarName = varName[i]
            longVarName  = varName[i]
            if longName != None: longVarName = longName
            var= rootgrp.createVariable(shortVarName,'f4',('time','lat','lon',) ,fill_value=MV,zlib=True)
            var.standard_name = varName[i]
            var.long_name = longVarName
            var.units = varUnits[i]
    else:    
        shortVarName = varName
        longVarName  = varName
        if longName != None: longVarName = longName
        var= rootgrp.createVariable(shortVarName,'f4',('time','lat','lon',) ,fill_value=MV,zlib=True)
        var.standard_name = varName
        var.long_name = longVarName
        var.units = varUnits
    rootgrp.title = title
    rootgrp.institution = "Radboud University, Utrecht University"
    rootgrp.processed = "Processed by Joyce Bosmans and Niko Wanders (joyce.bosmans@ru.nl and n.wanders@uu.nl)"
    rootgrp.sync()
    rootgrp.close()
    
def data2NetCDF(ncFile,varName,varField,timeStamp,posCnt = None):
  #-write data to netCDF
  rootgrp= nc.Dataset(ncFile,'a')    
  
  shortVarName= varName        
  
  date_time= rootgrp.variables['time']
  if posCnt == None: posCnt = len(date_time)
  
  date_time[posCnt]= nc.date2num(timeStamp,date_time.units,date_time.calendar)
  rootgrp.variables[shortVarName][posCnt,:,:]= (varField)
  
  rootgrp.sync()
  rootgrp.close()
  

def createNetCDFNoTime(ncFileName, varName, varUnits, latitudes, longitudes,\
                                      longName = None, loop=False, title="Bioclimoutput"):
    
    rootgrp= nc.Dataset(ncFileName,'w', format='NETCDF4')
    
    #-create dimensions - time is unlimited, others are fixed
    rootgrp.createDimension('lat',len(latitudes))
    rootgrp.createDimension('lon',len(longitudes))
    
    lat= rootgrp.createVariable('lat','f4',('lat',))
    lat.long_name= 'latitude'
    lat.units= 'degrees_north'
    lat.standard_name = 'latitude'
    
    lon= rootgrp.createVariable('lon','f4',('lon',))
    lon.standard_name= 'longitude'
    lon.long_name= 'longitude'
    lon.units= 'degrees_east'
    
    lat[:]= latitudes
    lon[:]= longitudes
    
    shortVarName = varName
    longVarName  = varName
    if longName != None: longVarName = longName
    var= rootgrp.createVariable(shortVarName,'f4',('lat','lon',) ,fill_value=MV,zlib=True)
    var.standard_name = varName
    var.long_name = longVarName
    var.units = varUnits
    rootgrp.title = title
    rootgrp.institution = "Radboud University, Utrecht University"
    rootgrp.processed = "Processed by Joyce Bosmans and Niko Wanders (joyce.bosmans@ru.nl and n.wanders@uu.nl)"
    
    rootgrp.sync()
    rootgrp.close()

def data2NetCDFNoTime(ncFile,varName,varField):
  #-write data to netCDF
  rootgrp= nc.Dataset(ncFile,'a')

  shortVarName= varName
  try:
    rootgrp.variables[shortVarName][:,:]= (varField)
  except:
    rootgrp.variables[shortVarName]= (varField)

  rootgrp.sync()
  rootgrp.close()
