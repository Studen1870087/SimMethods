#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 11:46:40 2019

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

"File Location"
path = '/Users/jake/Documents/CANES2019/Sim_Research/CANES_2019/files-tutorials/helix/runs_2-pulling/'
#path = '/Users/jake/Documents/CANES2019/Sim_Research/CANES_2019/files-tutorials/hairpin/runs_2_pulling/'
top = 'pull-v-0-01.tpr'
trj = 'pull-v-0-01.xtc'

"Load Universe"
u = mda.Universe(path+top, path+trj)

"Force-Distance data"
x_path = '/Users/jake/Documents/CANES2019/Sim_Research/Code/'
x_prot = 'helix'
x_file = '/force-vs-distance'
x_ext = '-v-0-01.dat'

"Load data"
data = np.loadtxt( x_path+x_prot+x_file+x_ext )
data[:,0] = -data[:,0]   #Helix only
lw=0.2

"""
"H Bonds from VMD Visualization (HAIRPIN)"
hydrogen_dict = {
    'H1' : u.select_atoms('bynum 140'),
    'H2' : u.select_atoms('bynum 153'),
    'H3' : u.select_atoms('bynum 33'),
    'H4' : u.select_atoms('bynum 123'),
    'H5' : u.select_atoms('bynum 63'),
    'H6' : u.select_atoms('bynum 89'),
    'H7' : u.select_atoms('bynum 50'),
}

nitrogen_dict = {
    'N1' : u.select_atoms('bonded bynum 140'),
    'N2' : u.select_atoms('bonded bynum 153'),
    'N3' : u.select_atoms('bonded bynum 33'),
    'N4' : u.select_atoms('bonded bynum 123'),
    'N5' : u.select_atoms('bonded bynum 63'),
    'N6' : u.select_atoms('bonded bynum 89'),
    'N7' : u.select_atoms('bonded bynum 50'),
}

oxygen_dict = {
    'O1' : u.select_atoms('bynum 164'),
    'O2' : u.select_atoms('bynum 10'),
    'O3' : u.select_atoms('bynum 130'),
    'O4' : u.select_atoms('bynum 40'),
    'O5' : u.select_atoms('bynum 100'),
    'O6' : u.select_atoms('bynum 71'),
    'O7' : u.select_atoms('bynum 76'),
}

oxygens = sum(list(oxygen_dict.values()))
hydrogens = sum(list(hydrogen_dict.values()))
nitrogens = sum(list(nitrogen_dict.values()))
"""

""
"H Bonds from VMD Visualization (HELIX)"
bond_dict = {
    'HB1' : u.select_atoms('bynum 137 or bynum 109'), #ASN8450:H & ASN8446:O
    'HB2' : u.select_atoms('bynum 128 or bynum 98'), #LEU8449:H & PHE8445:O
    'HB3' : u.select_atoms('bynum 119 or bynum 81'), #ASP8448:H & GLN8444:O
    'HB4' : u.select_atoms('bynum 111 or bynum 81'), #SER8447:H & GLN8444:O
    'HB5' : u.select_atoms('bynum 100 or bynum 57'), #ASN8446:H & TRP8442:O
    'HB6' : u.select_atoms('bynum 83 or bynum 36'), #PHE8445:H & LYS8441:O
    'HB7' : u.select_atoms('bynum 71 or bynum 23'), #GLN8444:H & GLN8440:O
    'HB8' : u.select_atoms('bynum 59 or bynum 23'), #GLN8443:H & GLN8440:O
    'HB9' : u.select_atoms('bynum 13 or bynum 64'), #GLN8440:H & GLN8443:OE1
    'HB10' : u.select_atoms('bynum 67 or bynum 18') #GLN8443:HE22 & GLN8440:OE1
}

hydrogen_dict = {
    'H1' : u.select_atoms('bynum 137'),
    'H2' : u.select_atoms('bynum 128'),
    'H3' : u.select_atoms('bynum 119'),
    'H4' : u.select_atoms('bynum 111'),
    'H5' : u.select_atoms('bynum 100'),
    'H6' : u.select_atoms('bynum 83'),
    'H7' : u.select_atoms('bynum 71'),
    'H8' : u.select_atoms('bynum 59'),
    'H9' : u.select_atoms('bynum 13'),
    'H10' : u.select_atoms('bynum 67')
}

nitrogen_dict = {
    'N1' : u.select_atoms('bonded bynum 137'),
    'N2' : u.select_atoms('bonded bynum 128'),
    'N3' : u.select_atoms('bonded bynum 119'),
    'N4' : u.select_atoms('bonded bynum 111'),
    'N5' : u.select_atoms('bonded bynum 100'),
    'N6' : u.select_atoms('bonded bynum 83'),
    'N7' : u.select_atoms('bonded bynum 71'),
    'N8' : u.select_atoms('bonded bynum 59'),
    'N9' : u.select_atoms('bonded bynum 13'),
    'N10' : u.select_atoms('bonded bynum 67')
}

oxygen_dict = {
    'O1' : u.select_atoms('bynum 109'),
    'O2' : u.select_atoms('bynum 98'),
    'O3' : u.select_atoms('bynum 81'),
    'O4' : u.select_atoms('bynum 81'),
    'O5' : u.select_atoms('bynum 57'),
    'O6' : u.select_atoms('bynum 36'),
    'O7' : u.select_atoms('bynum 23'),
    'O8' : u.select_atoms('bynum 23'),
    'O9' : u.select_atoms('bynum 64'),
    'O10' : u.select_atoms('bynum 18')
}

oxygens = sum(list(oxygen_dict.values()))
hydrogens = sum(list(hydrogen_dict.values()))
nitrogens = sum(list(nitrogen_dict.values()))
""

"Run HBond Autocorrelation"
hb_ac = HydrogenBondAutoCorrel(u, 
                               acceptors = oxygens,
                               hydrogens = hydrogens,
                               donors = nitrogens,
                               bond_type='continuous')
hb_ac.run()
hb_ac.solve()
tau = hb_ac.solution['tau']
time = hb_ac.solution['time']
results = hb_ac.solution['results']
estimate = hb_ac.solution['estimate']

"Max Range"
factor = 0.5
max_d = int(data[:,0].size * factor)
max_h = int(time.size * factor)
max_r = int(len(u.trajectory) * factor)

"Plot Results and Find Correlation"
fig, ax1 = plt.subplots()
color = 'red'
ax1.set_xlabel('Pull Distance ($\AA$)')
ax1.set_ylabel('Force (kJ/mol/nm)')
ax1.plot(data[0:max_d,0], data[0:max_d,1], '-.', lw=lw, color = color, label = 'Force')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twiny().twinx()

color = 'purple'
#ax2.set_xlabel('step')
ax2.set_xticks([])
ax2.set_ylabel('Autocorrelation')
ax2.plot(time[0:max_h], results[0:max_h], 'o', color=color)
ax2.plot(time[0:max_h], estimate[0:max_h], color=color, label = 'HBond Autocorr.')
ax2.tick_params(labelcolor=color)

xnew = np.linspace(0, max(time[0:max_h]), num=data[0:max_d,0].size, endpoint=True)
i_estimate = interp1d(time[0:max_h], estimate[0:max_h], 'slinear')

#ax2.plot(xnew, i_estimate(xnew), label='interpolated')

corr = pearsonr(data[0:max_d,1], i_estimate(xnew))

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

#ax1.axvline(x = data[3845,0], linewidth=4, color='b', ls =  '--')
#ax1.axvline(x = data[8520,0], linewidth=4, color='b', ls =  '--')

#plt.legend(lines1 + lines2, labels1 + labels2, loc=0, title='Correlation: '+"% 4.4f" %round(corr[0],4))
#fig.tight_layout()
fig.set_dpi(100)
plt.show()

"HBond Average Separation"
hbond_dst = []

for ts in u.trajectory:    
    tmp = 0
    for i in range(len(oxygens)):
        tmp += norm(oxygens[i].position - hydrogens[i].position)
    hbond_dst.append((u.trajectory.time, tmp/len(oxygens)))
    
hbond_dst = np.array(hbond_dst)


"Plot Results and Find Correlation"
fig, ax1 = plt.subplots()
color = 'red'
ax1.set_xlabel('Pull Distance ($\AA$)')
ax1.set_ylabel('Force (kJ/mol/nm)')
ax1.plot(data[0:max_d,0], data[0:max_d,1], '.', lw=lw, color = color, label = 'Force')
ax1.tick_params(axis='y', labelcolor=color)


ax2 = ax1.twiny().twinx()

color = 'purple'
#ax2.set_xlabel('step')
ax2.set_xticks([])
ax2.set_ylabel('N-H-O Distance ($\AA$)')
ax2.plot(hbond_dst[0:max_r,0], hbond_dst[0:max_r,1], color=color, label = 'Avg. N-H-O Separation')
ax2.tick_params(labelcolor=color)

xnew = np.linspace(0, max(hbond_dst[0:max_r,0]), num=data[0:max_d,0].size, endpoint=True)
i_hbond = interp1d(hbond_dst[0:max_r,0], hbond_dst[0:max_r,1], 'slinear')

#ax2.plot(xnew, i_hbond(xnew), label='interpolated')

corr = pearsonr(data[0:max_d,1], i_hbond(xnew))

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

#ax1.axvline(x = 1.8, linewidth=4, color='b', ls =  '--')
#ax1.axvline(x = 2.5, linewidth=4, color='b', ls =  '--')

plt.legend(lines1 + lines2, labels1 + labels2, loc=2, title='Correlation: '+"% 4.4f" %round(corr[0],4))
#fig.tight_layout()
fig.set_dpi(100)
plt.show()
