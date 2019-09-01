#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 15:53:50 2019

@author: jake
"""

import gromacs.formats
from matplotlib import gridspec
import matplotlib.pyplot as plt
import brewer2mpl

prot = 'hairpin'
ext = '.xvg'
data1 = gromacs.formats.XVG(prot+'/histo_2000'+ext)
data2 = gromacs.fileformats.xvg.XVG(filename=prot+'/bsResult_2000_trajG'+ext)
arr1 = data1.array
arr2 = data2.array

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

# FIGURE SETUP
fig = plt.figure(figsize=(11,8.5))
gs = gridspec.GridSpec(2,1,wspace=0.4,hspace=0.1)
ax1 = fig.add_subplot(gs[0,0])
clrs = [brewer2mpl.get_map('Paired','qualitative',9).mpl_colors[i] for i in range(9)]
clrs = clrs[4:] + clrs[:4] #Rotate colors

# PLOT HISTO
for i in range(len(data1.array)-1):
	ax1.plot(data1.array[0],zero_to_nan(data1.array[i+1]),c=clrs[i%8],lw=4)
	ax1.fill_between(data1.array[0],0,data1.array[i+1],facecolor=clrs[i%8],alpha=0.5)
#ax1.plot(data1.array[0], helix_ccnt, c=clrs[8], lw=5)
#ax1.plot(data1.array[0], hair_ccnt, c=clrs[8], lw=5)

ax1.grid()
ax1.set_ylabel('Counts')
ax1.get_yaxis().set_label_coords(-0.1,0.5)
#ax1.set_xlabel('Pull Distance ($\AA$)')
ax1.set_ylim([0,800])
#ax1.set_xticklabels([]) 
#ax1.set_ylim([0,-1.5])
#ax1.get_yaxis().set_label_coords(-0.25,0.5)
#ax1.yaxis.set_major_locator(MaxNLocator(5, prune='lower'))
#fig.set_dpi(100)

# PLOT PMF
ax2 = fig.add_subplot(gs[1,0])
#ax2 = ax1
ax2.plot(data2.array[0],data2.array[1],c=clrs[0],lw=4)
ax2.grid()
ax2.set_ylabel('Integrated PMF (kCal/mol)')
ax2.get_yaxis().set_label_coords(-0.1,0.5)
ax2.set_xlabel('Pull Distance ($\AA$)')
#ax1.set_xticklabels([]) 
#ax1.set_ylim([0,-1.5])
#ax1.get_yaxis().set_label_coords(-0.25,0.5)
#ax1.yaxis.set_major_locator(MaxNLocator(5, prune='lower'))

gs.tight_layout(fig, rect=[0.0, 0.0, 1, 1]) # Leave space for the common x-label.
#plt.close()
