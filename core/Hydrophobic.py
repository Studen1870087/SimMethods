#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:37:47 2019

@author: jake
"""

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda
from scipy.stats.stats import pearsonr
from scipy.interpolate import interp1d
from MDAnalysis.analysis.waterdynamics import HydrogenBondLifetimes as HBL
from MDAnalysis.analysis.waterdynamics import WaterOrientationalRelaxation as WOR
from MDAnalysis.analysis.waterdynamics import AngularDistribution as AD
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD
from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP


"File Location"
#path = '/Users/jake/Documents/CANES2019/Sim_Research/CANES_2019/files-tutorials/helix/runs_2-pulling/'
path = '/Users/jake/Documents/CANES2019/Sim_Research/CANES_2019/files-tutorials/hairpin/runs_2_pulling/'
top = 'pull.tpr'
trj = 'pull.xtc'

"Load Universe"
u = mda.Universe(path+top, path+trj)

"Force-Distance data"
x_path = '/Users/jake/Documents/CANES2019/Sim_Research/Code/'
x_prot = 'hairpin'
x_file = '/force-vs-distance'
x_ext = '.dat'

"Load data"
data = np.loadtxt( x_path+x_prot+x_file+x_ext )

""
"Water Dynamics"
u = mda.Universe(path+top, path+trj)
selection = "resname SOL and around 5.0 protein"

s1, e1 = 0, 10
s2, e2 = 323, 333

WOR_analysis1 = WOR(u, selection, s1, e1, e1-s1)
WOR_analysis1.run()
WOR_analysis2 = WOR(u, selection, s2, e2, e2-s2)
WOR_analysis2.run()

time = 0
for WOR_OH, WOR_HH, WOR_dip in WOR_analysis1.timeseries:
      print("{time} {WOR_OH} {time} {WOR_HH} {time} {WOR_dip}".format(time=time, WOR_OH=WOR_OH, WOR_HH=WOR_HH,WOR_dip=WOR_dip))
      time += 1
time = 0
for WOR_OH, WOR_HH, WOR_dip in WOR_analysis2.timeseries:
      print("{time} {WOR_OH} {time} {WOR_HH} {time} {WOR_dip}".format(time=time, WOR_OH=WOR_OH, WOR_HH=WOR_HH,WOR_dip=WOR_dip))
      time += 1

#now, if we want, we can plot our data
plt.figure(1,figsize=(18, 6))

#WOR OH
plt.subplot(131)
plt.xlabel('time')
plt.ylabel('WOR')
plt.title('WOR OH')
plt.plot(range(0,time),[column[0] for column in WOR_analysis1.timeseries])
plt.plot(range(0,time),[column[0] for column in WOR_analysis2.timeseries])

#WOR dip
plt.subplot(133)
plt.xlabel('time')
plt.ylabel('WOR')
plt.title('WOR dip')
plt.plot(range(0,time),[column[2] for column in WOR_analysis1.timeseries])
plt.plot(range(0,time),[column[2] for column in WOR_analysis2.timeseries])

plt.show()
""


"""
"Run Water Dynamics"
u = mda.Universe(path+top, path+trj)
selection1 = "resname SOL and around 5.0 protein"
selection2 = "protein"
s1 = 0
e1 = 10
s2 = len(u.trajectory)-10
e2 = len(u.trajectory)
HBL_analysis1 = HBL(u, selection1, selection2, s1, e1, e1-s1)
HBL_analysis1.run()
HBL_analysis2 = HBL(u, selection1, selection2, s2, e2, e2-s2)
HBL_analysis2.run()

time = 0
for HBLc, HBLi in HBL_analysis1.timeseries:
    print("{time} {HBLc} {time} {HBLi}".format(time=time, HBLc=HBLc, HBLi=HBLi))
    time += 1
time = 0
for HBLc, HBLi in HBL_analysis2.timeseries:
    print("{time} {HBLc} {time} {HBLi}".format(time=time, HBLc=HBLc, HBLi=HBLi))
    time += 1

#we can also plot our data
plt.figure(1,figsize=(18, 6))

#HBL continuos
plt.subplot(121)
plt.xlabel('time (ps)')
plt.ylabel('HBL')
plt.plot(range(0,time),[column[1] for column in HBL_analysis1.timeseries], label = 'Folded')
plt.plot(range(0,time),[column[1] for column in HBL_analysis2.timeseries], label = 'Unfolded')
plt.legend()

plt.show()
"""



"""
"Water Dynamics"
u = mda.Universe(path+top, path+trj)
selection = "resname SOL and around 5.0 protein"

s1, e1 = 0, 4
s2, e2 = 36, 40

MSD_analysis1 = MSD(u, selection, s1, e1, e1-s1)
MSD_analysis1.run()
MSD_analysis2 = MSD(u, selection, s2, e2, e2-s2)
MSD_analysis2.run()

time = 0
for msd in MSD_analysis1.timeseries:
      print("{time} {msd}".format(time=time, msd=msd))
      time += 1
time = 0
for msd in MSD_analysis2.timeseries:
      print("{time} {msd}".format(time=time, msd=msd))
      time += 1

#Plot
plt.xlabel('time')
plt.ylabel('MSD')
plt.title('MSD')
plt.plot(range(0,time),MSD_analysis1.timeseries)
plt.plot(range(0,time),MSD_analysis2.timeseries)
plt.show()
"""
