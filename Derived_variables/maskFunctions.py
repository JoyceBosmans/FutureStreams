import os
import netCDF4 as nc
import numpy as np
from functionsNC import *
import uuid
#import scipy.stats

### Created by Niko Wanders (n.wanders@uu.nl) and Joyce Bosmans (joyce.bosmans@ru.nl)

#os.system('module load cdo')	# not needed for masking?
indir   = '/vol/milkundata/FutureStreams'
outdir  = '/vol/milkundata/FutureStreams'

Qthreshold  = 10	# m3/s, grid cells with Q-mean < Q-threshold will be masked in Q-mask
WTthreshold = 350  	# K, grid cells with WT-threshold < WT-mean will be masked in WT-mask

########## General settings #############
futureStart   = [2021, 2041, 2061, 2081]
futureEnd     = [2040, 2060, 2080, 2099]

historicStart = [1976]
historicEnd   = [2005]
E2OStart      = [1979]

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
    fileName = '%s/weekAvg/%s_weekAvg_output_%s_%s_%d-01-07_to_%d-12-30.nc4' %(indir, varName,model,scen,start,end)
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
  if clim in ["Q-mask"]:
    varName = "discharge"
    for inputfile, start, end in zip(makeAnnualFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      inputfile  = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc4' %(outdir,'Q-mean','Q-mean',model,scen,start,end)
      outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31_%d.nc' %(outdir,clim,clim,model,scen,start,end,Qthreshold)
      print("%s %s %s %d" %(clim, model, scen, start))
      createMaskQ(inputfile,outputfile,Qthreshold)
  if clim in ["WT-mask"]:
    varName = "waterTemp"
    for inputfile, start, end in zip(makeAnnualFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      inputfile  = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc4' %(outdir,'WT-mean','WT-mean',model,scen,start,end)
      outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31_%d.nc' %(outdir,clim,clim,model,scen,start,end,WTthreshold)
      print("%s %s %s %d" %(clim, model, scen, start))
      createMaskWT(inputfile,outputfile,WTthreshold)


######################################################
def createMaskQ(inputfile,outputfile,threshold):
  data = getNCData(inputfile, varName = 'discharge')
  # set all land points to 1, then to nan if threshold is met. 
  mask = np.ones(data.shape)
  print(data.mean())
  mask = np.ma.masked_where(data < Qthreshold,mask)	#this already excludes ocean points
  latitudes  = getNCData(inputfile, varName = "latitude")
  longitudes = getNCData(inputfile, varName = "longitude")
  createNetCDFNoTime(outputfile, "mask", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile,"mask",mask)	# TO DO: check steps above, ocean now set to 1, land to nan?


def createMaskWT(inputfile,outputfile,threshold):
  data = getNCData(inputfile, varName = 'waterTemperature')
  print(type(data))
  print(data.shape)
  # set all land points to 1, then to nan if threshold is met. 
  mask = np.ones(data.shape)
  mask = np.ma.masked_where(data > WTthreshold,mask)
  
  latitudes  = getNCData(inputfile, varName = "latitude")
  longitudes = getNCData(inputfile, varName = "longitude")
  createNetCDFNoTime(outputfile, "mask", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile,"mask",mask)	# TO DO: check steps above, ocean now set to 1, land to nan?
