import os
import netCDF4 as nc
import numpy as np
from functionsNC import *
import uuid
import scipy.stats

### Created by Niko Wanders (n.wanders@uu.nl) and Joyce Bosmans (joyce.bosmans@ru.nl)

os.system('module load cdo')
indir   = '/scratch-shared/milkun'
outdir  = '/projects/0/milkun/output_derived/'

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

def makeMonthFileName(varName, model, scen):
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
    fileName = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, varName,model,scen,start,end)
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

#####################################################
######## Master function to call all data  #############################

def computeBioClimVariable(model, scen, clim):
  if clim in ["Q-mean"]:
    varName = "discharge"
    for inputfile, start, end in zip(makeAnnualFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,clim,clim,model,scen,start,end)
      print("%s %s %s %d" %(clim, model, scen, start))
      AnnualMean(inputfile,outputfile)
  if clim in ["WT-mean"]:
    varName = "waterTemp"
    for inputfile, start, end in zip(makeAnnualFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,clim,clim,model,scen,start,end)
      print("%s %s %s %d" %(clim, model, scen, start))
      AnnualMean(inputfile,outputfile)

  if clim in ['Q-max','Q-min','Q-range','Q-si','Q-wmin','Q-wmax']:
    varName = "discharge"
    for inputfile, start, end in zip(makeWeekFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      if clim == 'Q-min': 
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        Min(inputfile,outputfile)
      if clim == 'Q-max':
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        Max(inputfile,outputfile)
      if clim == 'Q-range':
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        inputfile1 = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, 'Q-max', 'Q-max',model,scen,start,end)
        inputfile2 = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, 'Q-min', 'Q-min',model,scen,start,end)
        AnnualRange(inputfile1,inputfile2, outputfile)
      if clim == 'Q-si':
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        inputfile2 = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,'Q-mean','Q-mean',model,scen,start,end)
        SeasonalityIndex(inputfile,inputfile2,outputfile)
      if clim == 'Q-wmin': 
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        TimeMin(inputfile,outputfile,varName)
      if clim == 'Q-wmax': 
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        TimeMax(inputfile,outputfile,varName)
  if clim in ['WT-max','WT-min','WT-range','WT-si','WT-wmax','WT-wmin']:
    varName = "waterTemp"
    for inputfile, start, end in zip(makeWeekFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      if clim == 'WT-min': 
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        Min(inputfile,outputfile)
      if clim == 'WT-max':
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        Max(inputfile,outputfile)
      if clim == 'WT-range':
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        inputfile1 = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, 'WT-max', 'WT-max',model,scen,start,end)
        inputfile2 = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, 'WT-min', 'WT-min',model,scen,start,end)
        AnnualRange(inputfile1,inputfile2, outputfile)
      if clim == 'WT-si':
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        inputfile2 = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,'WT-mean','WT-mean',model,scen,start,end)
        SeasonalityIndex(inputfile,inputfile2,outputfile)
      if clim == 'WT-wmin': 
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        TimeMin(inputfile,outputfile,varName)
      if clim == 'WT-wmax': 
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, clim, clim, model,scen,start,end)
        print("%s %s %s %d" %(clim, model, scen, start))
        TimeMax(inputfile,outputfile,varName)

  if clim in ["Q-zfw"]:
    varName = "discharge"
    for inputfile, start, end in zip(makeWeekFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      print("%s %s %s %d" %(clim, model, scen, start))
      outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,clim,clim,model,scen,start,end)
      ZeroFlow(inputfile,outputfile)
  if clim in ["WT-ztw"]:
    varName = "waterTemp"
    for inputfile, start, end in zip(makeWeekFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      print("%s %s %s %d" %(clim, model, scen, start))
      outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,clim,clim,model,scen,start,end)
      ZeroTemp(inputfile,outputfile)

  if clim in ["Q-bfi","Q-vi"]:
    varName = "discharge"
    for inputfile, start, end in zip(makeWeekFileName(varName, model, scen), getStartTimes(model, scen), getEndTimes(scen)):
      if clim == "Q-bfi":
        print("%s %s %s %d" %(clim, model, scen, start))
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,clim,clim,model,scen,start,end)
        baseFlow(inputfile,outputfile)
      if clim == "Q-vi":
        print("%s %s %s %d" %(clim, model, scen, start))
        outputfile = '%s/%s/%s_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,clim,clim,model,scen,start,end)
        variabilityIndex(inputfile,outputfile)



  if clim == "WettestDriestMonth":
    for start, end in zip(getStartTimes(model, scen), getEndTimes(scen)):
      print("%s %s %s %d" %(clim, model, scen, start))
      inputfile1  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'precipitation',model,scen,start,end)
      inputfile2  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'discharge',model,scen,start,end)
      inputfile3  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'waterTemp',model,scen,start,end)
      outputfile1 = '%s/BIOAUX/precipitation_wettest_month_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, model,scen,start,end)     # monthnumber of wettest month - AUX for auxilary files, shorten filename?
      outputfile2 = '%s/BIOAUX/precipitation_driest_month_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, model,scen,start,end)     # monthnumber of driest month - AUX for auxilary files, shorten filename?
      outputfile3 = '%s/Q-wm/Q-wm_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # discharge in wettest month
      outputfile4 = '%s/Q-dm/Q-dm_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)      # discharge in driest month
      outputfile5 = '%s/WT-wm/WT-wm_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)      # watertemp in wettest month
      outputfile6 = '%s/WT-dm/WT-dm_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)    # watertemp in driest month
      WettestDriestMonth(inputfile1,inputfile2,inputfile3,outputfile1,outputfile2,outputfile3,outputfile4,outputfile5,outputfile6)
  if clim == "HottestColdestMonth":
    for start, end in zip(getStartTimes(model, scen), getEndTimes(scen)):
      print("%s %s %s %d" %(clim, model, scen, start))
      inputfile1  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'tas',model,scen,start,end)
      inputfile2  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'discharge',model,scen,start,end)
      inputfile3  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'waterTemp',model,scen,start,end)
      outputfile1 = '%s/BIOAUX/tas_hottest_month_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, model,scen,start,end)
      outputfile2 = '%s/BIOAUX/tas_coldest_month_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, model,scen,start,end)
      outputfile3 = '%s/Q-hm/Q-hm_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # discharge in hottest month
      outputfile4 = '%s/Q-cm/Q-cm_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # discharge in coldest month
      outputfile5 = '%s/WT-hm/WT-hm_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # watertemp in hottest month
      outputfile6 = '%s/WT-cm/WT-cm_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # watertemp in coldest month
      HottestColdestMonth(inputfile1,inputfile2,inputfile3,outputfile1,outputfile2,outputfile3,outputfile4,outputfile5,outputfile6)
  if clim == "WettestDriestQuarter":
    for start, end in zip(getStartTimes(model, scen), getEndTimes(scen)):
      print("%s %s %s %d" %(clim, model, scen, start))
      inputfile1  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'precipitation',model,scen,start,end)
      inputfile2  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'discharge',model,scen,start,end)
      inputfile3  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'waterTemp',model,scen,start,end)
      outputfile1 = '%s/BIOAUX/precipitation_wettest_season_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, model,scen,start,end)
      outputfile2 = '%s/BIOAUX/precipitation_driest_season_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, model,scen,start,end)
      outputfile3 = '%s/Q-wq/Q-wq_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # discharge in hottest month
      outputfile4 = '%s/Q-dq/Q-dq_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # discharge in coldest month
      outputfile5 = '%s/WT-wq/WT-wq_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # watertemp in hottest month
      outputfile6 = '%s/WT-dq/WT-dq_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # watertemp in coldest month
      WettestDriestQuarter(inputfile1,inputfile2,inputfile3,outputfile1,outputfile2,outputfile3,outputfile4,outputfile5,outputfile6)
  if clim == "HottestColdestQuarter":
    for start, end in zip(getStartTimes(model, scen), getEndTimes(scen)):
      print("%s %s %s %d" %(clim, model, scen, start))
      inputfile1  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'tas',model,scen,start,end)
      inputfile2  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'discharge',model,scen,start,end)
      inputfile3  = '%s/%s_monthAvg_%s_%s_%d-01-31_to_%d-12-31.nc' %(indir, 'waterTemp',model,scen,start,end)
      outputfile1 = '%s/BIOAUX/tas_hottest_season_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, model,scen,start,end)
      outputfile2 = '%s/BIOAUX/tas_coldest_season_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir, model,scen,start,end)
      outputfile3 = '%s/Q-hq/Q-hq_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # discharge in hottest quarter
      outputfile4 = '%s/Q-cq/Q-cq_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # discharge in coldest quarter
      outputfile5 = '%s/WT-hq/WT-hq_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # watertemp in hottest quarter
      outputfile6 = '%s/WT-cq/WT-cq_%s_%s_%d-01-31_to_%d-12-31.nc' %(outdir,model,scen,start,end)     # watertemp in coldest quarter
      HottestColdestQuarter(inputfile1,inputfile2,inputfile3,outputfile1,outputfile2,outputfile3,outputfile4,outputfile5,outputfile6)

######################################################

def AnnualMean(inputfile,outputfile):
  os.system('cdo timmean %s %s' % (inputfile,outputfile))

def Min(inputfile,outputfile):
  randomFileName = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  os.system('cdo yearmin %s %s' % (inputfile, randomFileName))
  os.system('cdo timmean %s %s' % (randomFileName, outputfile))
  os.system('rm %s' %(randomFileName))

def Max(inputfile,outputfile):
  randomFileName = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  os.system('cdo yearmax %s %s' % (inputfile, randomFileName))
  os.system('cdo timmean %s %s' % (randomFileName, outputfile))
  os.system('rm %s' %(randomFileName))

def findMinMax(var1, var2):
  minMonth = np.argmin(var1, axis=0)
  minOutput = np.zeros((var1.shape[1], var2.shape[2]))
  maxMonth = np.argmax(var1, axis=0)
  maxOutput = np.zeros((var1.shape[1], var2.shape[2]))
  for i in range(var1.shape[0]):			# x.shape[0] => var1.shape[0]
    minSel = np.where(minMonth == i)
    minOutput[minSel] = var2[i,:,:][minSel]
    maxSel = np.where(maxMonth == i)
    maxOutput[maxSel] = var2[i,:,:][maxSel]		# maxOutput[minSel] => maxOutput[maxSel]
  return(minOutput, maxOutput)

def WettestDriestMonth(inputfile1,inputfile2,inputfile3,outputfile1,outputfile2,outputfile3,outputfile4,outputfile5,outputfile6):
  randomFileName1 = "temp/%s.nc.temp" %(str(uuid.uuid4()))	
  randomFileName2 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName3 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  os.system('cdo ymonmean %s %s' % (inputfile1, randomFileName1)) 	#precipitation
  os.system('cdo ymonmean %s %s' % (inputfile2, randomFileName2))       #discharge
  os.system('cdo ymonmean %s %s' % (inputfile3, randomFileName3))       #waterTemp
  data = getNCData(randomFileName1, varName = 'precipitation')
  data1 = getNCData(randomFileName2, varName = 'discharge')
  data2 = getNCData(randomFileName3, varName = 'waterTemperature')
  minData = np.argmin(data,axis=0)
  maxData = np.argmax(data,axis=0)
  minData1, maxData1 = findMinMax(data, data1)
  minData2, maxData2 = findMinMax(data, data2)
  latitudes = getNCData(randomFileName1, varName = "latitude")
  longitudes = getNCData(randomFileName1, varName = "longitude")
  createNetCDFNoTime(outputfile1, "monthNumber", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile1,"monthNumber",maxData)
  createNetCDFNoTime(outputfile2, "monthNumber", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile2,"monthNumber",minData)
  createNetCDFNoTime(outputfile3, "discharge", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile3,"discharge",maxData1)
  createNetCDFNoTime(outputfile4, "discharge", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile4,"discharge",minData1)
  createNetCDFNoTime(outputfile5, "waterTemperature", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile5,"waterTemperature",maxData2)
  createNetCDFNoTime(outputfile6, "waterTemperature", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile6,"waterTemperature",minData2)
  os.system('rm %s' %(randomFileName1))
  os.system('rm %s' %(randomFileName2))
  os.system('rm %s' %(randomFileName3))

def HottestColdestMonth(inputfile1,inputfile2,inputfile3,outputfile1,outputfile2,outputfile3,outputfile4,outputfile5,outputfile6):
  randomFileName1 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName2 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName3 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName4 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  gridfile = '/projects/0/milkun/output_derived/grid_5min.txt'
  os.system('cdo ymonmean %s %s' % (inputfile1, randomFileName1)) 	#airtemp
  os.system('cdo remapnn,%s %s %s' % (gridfile,randomFileName1,randomFileName2))
  os.system('cdo ymonmean %s %s' % (inputfile2, randomFileName3))       #discharge
  os.system('cdo ymonmean %s %s' % (inputfile3, randomFileName4))       #waterTemp
  if 'E2O' in inputfile1:
    data0 = getNCData(randomFileName2, varName = 'Tair')
    data  = data0[:,0,:,:]						#consider different varName and file structure
  else:
    data = getNCData(randomFileName2, varName = 'tas')
  data1 = getNCData(randomFileName3, varName = 'discharge')
  data2 = getNCData(randomFileName4, varName = 'waterTemperature')
  minData = np.argmin(data,axis=0)
  maxData = np.argmax(data,axis=0)
  minData1, maxData1 = findMinMax(data, data1)
  minData2, maxData2 = findMinMax(data, data2)
  latitudes = getNCData(randomFileName3, varName = "latitude")
  longitudes = getNCData(randomFileName3, varName = "longitude")
  createNetCDFNoTime(outputfile1, "monthNumber", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile1,"monthNumber",maxData)
  createNetCDFNoTime(outputfile2, "monthNumber", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile2,"monthNumber",minData)
  createNetCDFNoTime(outputfile3, "discharge", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile3,"discharge",maxData1)
  createNetCDFNoTime(outputfile4, "discharge", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile4,"discharge",minData1)
  createNetCDFNoTime(outputfile5, "waterTemperature", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile5,"waterTemperature",maxData2)
  createNetCDFNoTime(outputfile6, "waterTemperature", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile6,"waterTemperature",minData2)
  os.system('rm %s' %(randomFileName1))
  os.system('rm %s' %(randomFileName2))
  os.system('rm %s' %(randomFileName3))
  os.system('rm %s' %(randomFileName4))

def WettestDriestQuarter(inputfile1,inputfile2,inputfile3,outputfile1,outputfile2,outputfile3,outputfile4,outputfile5,outputfile6):
  randomFileName1 = "temp/%s.nc.temp" %(str(uuid.uuid4()))	
  randomFileName2 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName3 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  os.system('cdo yseasmean %s %s' % (inputfile1, randomFileName1)) 	#precipitation
  os.system('cdo yseasmean %s %s' % (inputfile2, randomFileName2))       #discharge
  os.system('cdo yseasmean %s %s' % (inputfile3, randomFileName3))       #waterTemp
  data = getNCData(randomFileName1, varName = 'precipitation')
  data1 = getNCData(randomFileName2, varName = 'discharge')
  data2 = getNCData(randomFileName3, varName = 'waterTemperature')
  minData = np.argmin(data,axis=0)
  maxData = np.argmax(data,axis=0)
  minData1, maxData1 = findMinMax(data, data1)
  minData2, maxData2 = findMinMax(data, data2)
  latitudes = getNCData(randomFileName1, varName = "latitude")
  longitudes = getNCData(randomFileName1, varName = "longitude")
  createNetCDFNoTime(outputfile1, "seasonNumber", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile1,"seasonNumber",maxData)
  createNetCDFNoTime(outputfile2, "seasonNumber", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile2,"seasonNumber",minData)
  createNetCDFNoTime(outputfile3, "discharge", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile3,"discharge",maxData1)
  createNetCDFNoTime(outputfile4, "discharge", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile4,"discharge",minData1)
  createNetCDFNoTime(outputfile5, "waterTemperature", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile5,"waterTemperature",maxData2)
  createNetCDFNoTime(outputfile6, "waterTemperature", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile6,"waterTemperature",minData2)
  os.system('rm %s' %(randomFileName1))
  os.system('rm %s' %(randomFileName2))
  os.system('rm %s' %(randomFileName3))

def HottestColdestQuarter(inputfile1,inputfile2,inputfile3,outputfile1,outputfile2,outputfile3,outputfile4,outputfile5,outputfile6):
  randomFileName1 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName2 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName3 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName4 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  gridfile = '/projects/0/milkun/output_derived/grid_5min.txt'
  os.system('cdo yseasmean %s %s' % (inputfile1, randomFileName1)) 	#airtemp
  os.system('cdo remapnn,%s %s %s' % (gridfile,randomFileName1,randomFileName2))
  os.system('cdo yseasmean %s %s' % (inputfile2, randomFileName3))       #discharge
  os.system('cdo yseasmean %s %s' % (inputfile3, randomFileName4))       #waterTemp
  if 'E2O' in inputfile1:
    data0 = getNCData(randomFileName2, varName = 'Tair')
    print('E2O Tair:', data0.shape)
    data  = data0[:,0,:,:]
    print('E2O Tair:', data.shape)
  else:
    data = getNCData(randomFileName2, varName = 'tas')
  data1 = getNCData(randomFileName3, varName = 'discharge')
  data2 = getNCData(randomFileName4, varName = 'waterTemperature')
  minData = np.argmin(data,axis=0)
  maxData = np.argmax(data,axis=0)
  minData1, maxData1 = findMinMax(data, data1)
  minData2, maxData2 = findMinMax(data, data2)
  latitudes = getNCData(randomFileName3, varName = "latitude")
  longitudes = getNCData(randomFileName3, varName = "longitude")
  createNetCDFNoTime(outputfile1, "seasonNumber", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile1,"seasonNumber",maxData)
  createNetCDFNoTime(outputfile2, "seasonNumber", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile2,"seasonNumber",minData)
  createNetCDFNoTime(outputfile3, "discharge", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile3,"discharge",maxData1)
  createNetCDFNoTime(outputfile4, "discharge", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile4,"discharge",minData1)
  createNetCDFNoTime(outputfile5, "waterTemperature", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile5,"waterTemperature",maxData2)
  createNetCDFNoTime(outputfile6, "waterTemperature", "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile6,"waterTemperature",minData2)
  os.system('rm %s' %(randomFileName1))
  os.system('rm %s' %(randomFileName2))
  os.system('rm %s' %(randomFileName3))
  os.system('rm %s' %(randomFileName4))

def ZeroFlow(inputfile,outputfile):
  randomFileName1 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName2 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  if 'hist'in inputfile:
    os.system('cdo lec,0.001 %s %s' % (inputfile, randomFileName1))
    os.system('cdo timsum %s %s' %(randomFileName1, randomFileName2))
    os.system('cdo divc,30 %s %s' %(randomFileName2, outputfile))
  if 'rcp'in inputfile:
    div = 20
    if '2081' in inputfile: div = 19
    os.system('cdo lec,0.001 %s %s' % (inputfile, randomFileName1))
    os.system('cdo timsum %s %s' %(randomFileName1, randomFileName2))
    os.system('cdo divc,%d %s %s' %(div, randomFileName2, outputfile))
  os.system('rm %s' %(randomFileName1))
  os.system('rm %s' %(randomFileName2))

def ZeroTemp(inputfile,outputfile):
  randomFileName1 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  randomFileName2 = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  if 'hist'in inputfile:
    os.system('cdo lec,273.55 %s %s' % (inputfile, randomFileName1))
    os.system('cdo timsum %s %s' %(randomFileName1, randomFileName2))
    os.system('cdo divc,30 %s %s' %(randomFileName2, outputfile))
  if 'rcp'in inputfile:
    div = 20
    if '2081' in inputfile: div = 19
    os.system('cdo lec,273.55 %s %s' % (inputfile, randomFileName1))
    os.system('cdo timsum %s %s' %(randomFileName1, randomFileName2))
    os.system('cdo divc,%d %s %s' %(div, randomFileName2, outputfile))
  os.system('rm %s' %(randomFileName1))
  os.system('rm %s' %(randomFileName2))

def AnnualRange(inputfile1,inputfile2, outputfile):
  if os.path.exists(inputfile1) and os.path.exists(inputfile2):
    os.system('cdo sub %s %s %s' % (inputfile1, inputfile2, outputfile))
  else:
    print('Either %s and or %s not available to compute AnnualRange'%(inputfile1,inputfile2))

def SeasonalityIndex(inputfile1,inputfile2,outputfile):
  randomFileName = "temp/%s.nc.temp" %(str(uuid.uuid4()))
  if os.path.exists(inputfile2):
    os.system('cdo timstd %s %s' % (inputfile1,randomFileName))
    os.system('cdo div %s %s %s' % (randomFileName, inputfile2, outputfile))
  else:
    print('%s not available to compute SeasonalityIndex'%(inputfile1))
  os.system('rm %s' %(randomFileName))

def TimeMin(inputfile,outputfile,varName):
  length = 52
  times = getNCData(inputfile, varName = 'time')
  timeValues = readTimeValues(inputfile)
  nYears = int(len(times)/length)
  minTime = np.zeros((nYears, 2160, 4320))
  if varName == 'waterTemp':
    varName = 'waterTemperature'
  for start, end, count in zip(range(0,len(times),52), range(51, len(times), 52), range(nYears)):
    data = getNCIndexData(inputfile, varName = varName, idxStart=start, idxEnd=end+1)
    minTime[count] = np.argmin(data, axis=0)
  modeData = scipy.stats.mode(minTime)
  sel = (modeData[1] == 1)
  output = modeData[0]
  output[sel] = np.nan
  latitudes = getNCData(inputfile, varName = "latitude")
  longitudes = getNCData(inputfile, varName = "longitude")
  createNetCDFNoTime(outputfile, varName, "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile,varName,output)

def TimeMax(inputfile,outputfile,varName):
  length = 52
  times = getNCData(inputfile, varName = 'time')
  timeValues = readTimeValues(inputfile)
  nYears = int(len(times)/length)
  maxTime = np.zeros((nYears, 2160, 4320))
  if varName == 'waterTemp':
    varName = 'waterTemperature'
  for start, end, count in zip(range(0,len(times),52), range(51, len(times), 52), range(nYears)):
    data = getNCIndexData(inputfile, varName = varName, idxStart=start, idxEnd=end+1)
    maxTime[count] = np.argmax(data, axis=0)
  modeData = scipy.stats.mode(maxTime)
  sel = (modeData[1] == 1)
  output = modeData[0]
  output[sel] = np.nan
  latitudes = getNCData(inputfile, varName = "latitude")
  longitudes = getNCData(inputfile, varName = "longitude")
  createNetCDFNoTime(outputfile, varName, "-", latitudes, longitudes)
  data2NetCDFNoTime(outputfile,varName,output)

def baseFlow(inputfile, outputfile):
  times = getNCData(inputfile, varName = 'time')
  timeValues = readTimeValues(inputfile)
  nYears = int(len(times)/52)
  outputQ90 = np.zeros((nYears, 2160, 4320))
  outputMean = np.zeros((nYears, 2160, 4320))
  for start, end, count in zip(range(0,len(times),52), range(51, len(times), 52), range(nYears)):
    allData = getNCIndexData(inputfile, varName = 'discharge', idxStart=start, idxEnd=end+1)
    outputQ90[count,:,:] = np.percentile(allData, 10, axis=0)
    outputMean[count,:,:] = np.mean(allData, axis=0)
  outputMean = np.where(outputMean == 0, np.nan, outputMean)	#added to avoid div by 0
  output = outputQ90/outputMean
  latitudes = getNCData(inputfile, varName = "latitude")
  longitudes = getNCData(inputfile, varName = "longitude")
  createNetCDF(outputfile, "baseFlowIndex", "-", latitudes, longitudes)
  for start, count in zip(range(0,len(times),52), range(nYears)):
    data2NetCDF(outputfile,"baseFlowIndex",output[count,:,:],timeValues[start])
  
def variabilityIndex(inputfile, outputfile):
  times = getNCData(inputfile, varName = 'time')
  timeValues = readTimeValues(inputfile)
  nYears = int(len(times)/52)
  outputQ25 = np.zeros((nYears, 2160, 4320))
  outputQ50 = np.zeros((nYears, 2160, 4320))
  outputQ75 = np.zeros((nYears, 2160, 4320))
  for start, end, count in zip(range(0,len(times),52), range(51, len(times), 52), range(nYears)):
    allData = getNCIndexData(inputfile, varName = 'discharge', idxStart=start, idxEnd=end+1)
    outputQ75[count,:,:] = np.percentile(allData, 25, axis=0)	# 75 -> 25
    outputQ25[count,:,:] = np.percentile(allData, 75, axis=0)	# 25 -> 72
    outputQ50[count,:,:] = np.percentile(allData, 50, axis=0)
  outputQ50 = np.where(outputQ50 == 0, np.nan, outputQ50)	#added to avoid div by 0
  output = (outputQ25-outputQ75)/outputQ50
  latitudes = getNCData(inputfile, varName = "latitude")
  longitudes = getNCData(inputfile, varName = "longitude")
  createNetCDF(outputfile, "variabilityIndex", "-", latitudes, longitudes)
  for start, count in zip(range(0,len(times),52), range(nYears)):
    data2NetCDF(outputfile,"variabilityIndex",output[count,:,:],timeValues[start])

