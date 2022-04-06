import os
import netCDF4 as nc
import numpy as np
from functionsNC import *
import uuid
#import scipy.stats

### Created by Niko Wanders (n.wanders@uu.nl) and Joyce Bosmans (joyce.bosmans@ru.nl)

#os.system('module load cdo')	# not needed for masking?
indir   = '/vol/milkundata/FutureStreams'		# check function makeWeekFileName for filenames including directory
outdir  = '/vol/milkundata/FutureStreams'		# check lines 82, 88 for filenames including directory

Qthreshold  = 10	# m3/s, grid cells with Q-mean < Q-threshold will be masked in Q-mask
WTthreshold = 350  	# K, grid cells with WT-threshold < WT-mean will be masked in WT-mask

########## General settings #############
futureStart   = [2006, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]
futureEnd     = [2019, 2029, 2039, 2049, 2059, 2069, 2079, 2089, 2099]

historicStart = [1976, 1986, 1996]
historicEnd   = [1985, 1995, 2005]
E2OStart      = [1979, 1986, 1996]

########### General functions ############
def makeAnnualFileName(varName, model, scen):
  if scen == "hist" and model == "E2O":
    startS = E2OStart
    endS = historicEnd
  elif scen == "hist":
    startS = historicStart
    endS = historicEnd
  else:
    startS = futureStart
    endS = futureEnd
  fileNames = []
  for start, end in zip(startS, endS):
    if varName == "discharge":
      fileName = '%s/%s_annuaAvg_%s_%s_%d-01-07_to_%d-12-30.nc' %(indir, varName,model,scen,start,end)
    elif varName == "waterTemp":
      fileName = '%s/%s_annuaAvg_%s_%s_%d-01-07_to_%d-12-30.nc' %(indir, varName,model,scen,start,end)
    fileNames.append(fileName)
  return(fileNames)

def makeWeekFileName(varName, model, scen):
  if scen == "hist" and model == "E2O":
    startS = E2OStart
    endS = historicEnd
  elif scen == "hist":
    startS = historicStart
    endS = historicEnd
  else:
    startS = futureStart
    endS = futureEnd
  fileNames = []
  for start, end in zip(startS, endS):
    fileName = '%s/weekAvg/%s_weekAvg_output_%s_%s_%d-01-07_to_%d-12-30.nc' %(indir, varName,model,scen,start,end)
    fileNames.append(fileName)
  return(fileNames)
  
def getStartTimes(model, scen):
  if scen == "hist" and model == "E2O":
    startS = E2OStart
  elif scen == "hist":
    startS = historicStart
  else:
    startS = futureStart
  return(startS)

def getEndTimes(scen):
  if scen == "hist":
    endS = historicEnd
  else:
    endS = futureEnd
  return(endS)

######## Master function to call all data  #############################
def computeBioClimVariable(model, scen, clim):
  if clim in ["Q-mask-weekly"]:
    varName = "discharge"
    for inputfile, start, end in zip(makeWeekFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      outputfile = '%s/%s/%s_weekAvg_%s_%s_%d-01-31_to_%d-12-31_masked_%d.nc' %(outdir,clim,varName,model,scen,start,end,Qthreshold)
      print("%s %s %s %d" %(clim, model, scen, start))
      createMaskQ(inputfile,outputfile,Qthreshold)
  if clim in ["WT-mask-weekly"]:
    varName = "waterTemp"
    for inputfile, start, end in zip(makeWeekFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      outputfile = '%s/%s/%s_weekAvg_%s_%s_%d-01-31_to_%d-12-31_masked_%d.nc' %(outdir,clim,varName,model,scen,start,end,WTthreshold)
      print("%s %s %s %d" %(clim, model, scen, start))
      createMaskWT(inputfile,outputfile,WTthreshold)


######################################################
def createMaskQ(inputfile,outputfile,threshold):
  data = getNCData(inputfile, varName = 'discharge')
  # for weekly files, directly apply mask to data instead of storing mask separately
  data = np.ma.masked_where(data < Qthreshold,data)	#this already excludes ocean points
  latitudes  = getNCData(inputfile, varName = "latitude")
  longitudes = getNCData(inputfile, varName = "longitude")
  timeValues = readTimeValues(inputfile)
  createNetCDF(outputfile, "discharge_masked", "m3.s-1", latitudes, longitudes)
  for i in range(0,len(timeValues)):
    print(timeValues[i])
    data2NetCDF(outputfile,"discharge_masked",data[i,:,:],timeValues[i])


def createMaskWT(inputfile,outputfile,threshold):
  data = getNCData(inputfile, varName = 'waterTemperature')	# in K
  # for weekly files, directly apply mask to data instead of storing mask separately
  data = np.ma.masked_where(data > WTthreshold,data)	#this already excludes ocean points  
  latitudes  = getNCData(inputfile, varName = "latitude")
  longitudes = getNCData(inputfile, varName = "longitude")
  timeValues = readTimeValues(inputfile)
  createNetCDF(outputfile, "waterTemperature_masked", "K", latitudes, longitudes)
  for i in range(0,len(timeValues)):
    #print(i,timeValues[i])
    data2NetCDF(outputfile,"waterTemperature_masked",data[i,:,:],timeValues[i])
