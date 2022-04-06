import os,sys
import netCDF4 as nc
import numpy as np
import maskFunctionsWeekly as mask
from functionsNC import *
import multiprocessing as mp

### Created by Niko Wanders (n.wanders@uu.nl) and Joyce Bosmans (joyce.bosmans@ru.nl)
### in- and output directories and thresholds are set in maskFunctions.py. Subdirectories per variable in the output directory are to be created prior to running these scripts. 

modelS        = ["gfdl", "hadgem", "ipsl", "miroc", "noresm","E2O"]	
scenS         = ["hist", "rcp2p6", "rcp4p5", "rcp6p0", "rcp8p5"]

# ~ ###alternatively, specify model, scenario in command line:
# ~ modelS        = [sys.argv[1]]
# ~ scenS         = [sys.argv[2]]

climVars = ['Q-mask-weekly','WT-mask-weekly']

mS = []
sS = []
cS = []

for model in modelS:
  if model == 'E2O':
    scen = 'hist'
    for clim in climVars:
      mS.append(model)
      sS.append(scen)
      cS.append(clim)

  else:
    for scen in scenS:
      for clim in climVars:
        mS.append(model)
        sS.append(scen)
        cS.append(clim)

### lines below for testing at login node / without parallel processes:
for num in range(len(cS)):
  mask.computeBioClimVariable(mS[num], sS[num], cS[num])


### lines below for parallel processes:
#pool = mp.Pool(processes=2)

#results2 = [pool.apply_async(bio.computeBioClimVariable,args=(c, m, s)) for c,m,s in zip(cS, mS, sS)]
#outputList2 = [p.get() for p in results2]



