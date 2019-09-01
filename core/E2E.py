#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 20:50:52 2019

@author: jake
"""
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.hbonds import HydrogenBondAnalysis, HydrogenBondAutoCorrel
import MDAnalysis.analysis.rdf as MDRDF
from scipy.stats.stats import pearsonr
from scipy.interpolate import interp1d
from scipy.stats import norm as norms
import matplotlib.mlab as mlab

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

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
data[:,0] = -data[:,0]
lw=0.2

"Select Atoms"
putin = u.select_atoms('protein')
endCA = u.select_atoms('(resid 0 or resid 11) and name CA') #CHANGE RESID FOR HAIRPIN/HELIX

Rgyr = []
e2e = []
for ts in u.trajectory:
   Rgyr.append((u.trajectory.time, putin.radius_of_gyration()))
   e2e.append((u.trajectory.time, norm(endCA[0].position - endCA[1].position)))
e2e = np.array(e2e)
Rgyr = np.array(Rgyr)

"PLOT THAT SHIT"
fig, ax1 = plt.subplots(dpi=100)
color = 'green'
ax1.set_xlabel('time (ps)')
ax1.set_ylabel('Radius of Gyration ($\AA$)')
ax1.plot(Rgyr[:,0], Rgyr[:,1], color = color, label = 'RoG')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()

color = 'purple'
#ax2.set_xlabel('step')
ax2.set_ylabel('End-to-end distance ($\AA$)')
ax2.plot(e2e[:,0], e2e[:,1], ':', color = color, label = 'E-to-E dist.')
ax2.tick_params(labelcolor=color)

xnew = np.linspace(0, max(Rgyr[:,0]), num=data[:,0].size, endpoint=True)
i_estimate = interp1d(Rgyr[:,0], Rgyr[:,1], 'slinear')

corr = pearsonr(i_estimate(xnew), data[:,1])

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

plt.legend(lines1 + lines2, labels1 + labels2, loc=2,
           title='Correlation w/ Force: '+"% 4.4f" %round(corr[0],4))


rect = [0.5,0.14,0.35,0.35]
ax3 = add_subplot_axes(ax2,rect)

color = 'green'
ax3.set_xlabel('time (ps)')
ax3.set_ylabel('RoG ($\AA$)')
ax3.plot(Rgyr[:,0], Rgyr[:,1], color = color, label = 'RoG')
ax3.tick_params(axis='y', labelcolor=color)

ax4 = ax3.twinx()

color = 'red'
newx = np.linspace(0, max(Rgyr[:,0]), len(data[:,1]))
ax4.plot(newx, data[:,1], '-.', color = color, label = 'Force')
ax4.set_ylabel('Force (kJ/mol/nm)')
ax4.tick_params(axis='y', labelcolor=color)

#fig.tight_layout()
plt.show()
"""

"File Location"
path = '/Users/jake/Documents/CANES2019/Sim_Research/CANES_2019/files-tutorials/helix/runs_2-pulling/'
#path = '/Users/jake/Documents/CANES2019/Sim_Research/CANES_2019/files-tutorials/hairpin/runs_2_pulling/'
top1 = 'pull-v-0-1.tpr'
trj1 = 'pull-v-0-1.xtc'
top2 = 'pull-v-0-01.tpr'
trj2 = 'pull-v-0-01.xtc'
top3 = 'pull-v-0-001.tpr'
trj3 = 'pull-v-0-001.xtc'
top4 = 'pull-v-0-0001.tpr'
trj4 = 'pull-v-0-0001.xtc'

"Load Universe"
u1 = mda.Universe(path+top1, path+trj1)
u2 = mda.Universe(path+top2, path+trj2)
u3 = mda.Universe(path+top3, path+trj3)
u4 = mda.Universe(path+top4, path+trj4)

"Select Atoms"
endCA1 = u1.select_atoms('(resid 0 or resid 12) and name CA') #CHANGE RESID FOR HAIRPIN/HELIX
endCA2 = u2.select_atoms('(resid 0 or resid 12) and name CA') #CHANGE RESID FOR HAIRPIN/HELIX
endCA3 = u3.select_atoms('(resid 0 or resid 12) and name CA') #CHANGE RESID FOR HAIRPIN/HELIX
endCA4 = u4.select_atoms('(resid 0 or resid 12) and name CA') #CHANGE RESID FOR HAIRPIN/HELIX


e2e1 = []
e2e2 = []
e2e3 = []
e2e4 = []
for ts in u1.trajectory:
    e2e1.append((u1.trajectory.time, norm(endCA1[0].position - endCA1[1].position)))
for ts in u2.trajectory:
    e2e2.append((u2.trajectory.time, norm(endCA2[0].position - endCA2[1].position)))
for ts in u3.trajectory:
    e2e3.append((u3.trajectory.time, norm(endCA3[0].position - endCA3[1].position)))
for ts in u4.trajectory:
    e2e4.append((u4.trajectory.time, norm(endCA4[0].position - endCA4[1].position)))

e2e1 = np.array(e2e1)
e2e2 = np.array(e2e2)
e2e3 = np.array(e2e3)
e2e4 = np.array(e2e4)

fig = plt.subplots(dpi=100)

xran = np.linspace(10,60,1000)

mu, sigma = norms.fit(e2e1[:,1])
n, bins, patches = plt.hist(e2e1[:,1], normed=1, bins = 10, label = '0.1 nm/ps', alpha = 0.3)
y = norms.pdf( xran, mu, sigma)
l = plt.plot(xran, y, 'b--', linewidth=2)

(mu, sigma) = norms.fit(e2e2[:,1])
n, bins, patches = plt.hist(e2e2[:,1], normed=1, bins = 10, label = '0.01 nm/ps', color = 'm', alpha = 0.3)
y = norms.pdf( xran, mu, sigma)
l = plt.plot(xran, y, 'm--', linewidth=2)

(mu, sigma) = norms.fit(e2e3[:,1])
n, bins, patches = plt.hist(e2e3[:,1], normed=1, bins = 10, label = '0.001 nm/ps', color = 'g', alpha = 0.3)
y = norms.pdf( xran, mu, sigma)
l = plt.plot(xran, y, 'g--', linewidth=2)

(mu, sigma) = norms.fit(e2e4[:,1])
n, bins, patches = plt.hist(e2e4[:,1], normed=1, bins = 10, label = '0.0001 nm/ps', color = 'r', alpha = 0.3)
y = norms.pdf( xran, mu, sigma)
l = plt.plot(xran, y, 'r--', linewidth=2)

plt.xlabel('End-to-end distance $\AA$')
plt.ylabel('Frequency')
plt.grid()

plt.legend()