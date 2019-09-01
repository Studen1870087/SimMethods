#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 18:10:46 2019

@author: jake
"""

import gromacs.formats
from matplotlib import gridspec
import matplotlib.pyplot as plt
import brewer2mpl

prot = 'hairpin'
ext = '_bhist.xvg'

data1 = gromacs.formats.XVG(prot+'/bsResult_2000'+ext)
arr1 = data1.array

# FIGURE SETUP
fig = plt.figure(figsize=(11,8.5))
gs = gridspec.GridSpec(2,1,wspace=0.4,hspace=0.1)
ax1 = fig.add_subplot(gs[0,0])
clrs = [brewer2mpl.get_map('Paired','qualitative',9).mpl_colors[i] for i in range(9)]
clrs = clrs[4:] + clrs[:4] #Rotate colors

# PLOTTING
ax2 = ax1
ax2.errorbar(data1.array[0],data1.array[1], data1.array[2],c=clrs[0],lw=1)
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