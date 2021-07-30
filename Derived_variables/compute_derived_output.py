import os,sys
import netCDF4 as nc
import numpy as np
import derivedFunctions as bio
from functionsNC import *
import multiprocessing as mp

### Created by Niko Wanders (n.wanders@uu.nl) and Joyce Bosmans (joyce.bosmans@ru.nl)
### in- and output directories are set in derivedFunctions.py. Subdirectories per variable in the output directory are to be set prior to running these scripts. 

os.system('module load cdo')

modelS        = ["gfdl", "hadgem", "ipsl", "miroc", "noresm","E2O"]	
scenS         = ["hist", "rcp2p6", "rcp4p5", "rcp6p0", "rcp8p5"]

###alternatively, specify model in command line:
#modelS        = [sys.argv[1]]

climVars = ['Q-mean','Q-min','Q-max','WettestDriestMonth','HottestColdestMonth','WettestDriestQuarter','HottestColdestQuarter','Q-zfw','Q-range','Q-si','Q-bfi','Q-hvi','Q-wmin','Q-wmax','WT-mean','WT-min','WT-max','WT-ztw','WT-range','WT-si','WT-wmin','WT-wmax']
### total: 16 WT and 18 QT variables (see Tables 2 and 3), some are grouped:
### WettestDriestMonth:    Q-wm, Q-dm, WT-wm, WT-dm
### HottestColdestMonth:   Q-hm, Q-cm, WT-hm, WT-cm
### WettestDriestQuarter:  Q-wq, Q-dq, WT-wq, WT-dq
### HottestColdestQuarter: Q-hq, Q-cq, WT-hq, WT-cq

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
  bio.computeBioClimVariable(mS[num], sS[num], cS[num])


### lines below for parallel processes:
#pool = mp.Pool(processes=2)

#results2 = [pool.apply_async(bio.computeBioClimVariable,args=(c, m, s)) for c,m,s in zip(cS, mS, sS)]
#outputList2 = [p.get() for p in results2]


