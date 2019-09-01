 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 12:31:47 2019

@author: jake
"""

import numpy as np
import gromacs.formats
import matplotlib.pyplot as plt
from heatmap import corrplot
import pandas as pd
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 100

hair_path = 'hairpin/'
hair_ext = '_traj.xvg'
helix_path = 'helix/'
helix_ext = '_ac.xvg'

hair_d = gromacs.formats.XVG(hair_path+'histo_20'+hair_ext)
helix_d = gromacs.formats.XVG(helix_path+'histo_20'+helix_ext)

hair_arr = hair_d.array.T
helix_arr = helix_d.array.T

def histogram_intersection(h1, h2):
    sm = 0
    h1 /= sum(h1)
    h2 /= sum(h2)
    for i in range(len(h1)):
        sm += min(h1[i], h2[i])
    return sm

maxsort = np.array([])
for i in range(1, len(helix_d.array)):
    maxsort = np.append(maxsort, np.argmax(helix_d.array[i]))

new_helix = np.zeros((len(helix_d.array), len(helix_d.array[1])))
for i in range(1,8):
    new_helix[i] = helix_d.array[i]
new_helix[i+1] = helix_d.array[10]
new_helix[i+2] = helix_d.array[8]
new_helix[i+3] = helix_d.array[9]
new_helix[i+4] = helix_d.array[11]
new_helix[i+5] = helix_d.array[13]
new_helix[i+6] = helix_d.array[12]
for i in range(14,21):
    new_helix[i] = helix_d.array[i]
new_helix[i+1] = helix_d.array[22]
new_helix[i+2] = helix_d.array[21]
new_helix[i+3] = helix_d.array[23]
new_helix[i+4] = helix_d.array[25]
new_helix[i+5] = helix_d.array[24]

new_helix = new_helix.T
new_helix = np.delete(new_helix, 0, axis=1)

maxsort = np.array([])
for i in range(1, len(new_helix)):
    maxsort = np.append(maxsort, np.argmax(new_helix[i]))


mat = np.zeros((len(hair_d.array)-1, len(hair_d.array)-1))
for i in range(1, len(hair_d.array)):
    for j in range(1,i):
        mat[i-1,j-1] = histogram_intersection(hair_d.array[i], hair_d.array[j])
        mat[j-1,i-1] = mat[i-1,j-1]
        mat[i-1,i-1] = 1
mat[0,0] = 1
plt.imshow(mat)
plt.show()

df = pd.DataFrame(mat)
corrplot(df)
plt.show()




mat = np.zeros((len(new_helix.T), len(new_helix.T)))
for i in range(0, len(new_helix.T)):
    for j in range(0,i):
        mat[i,j] = histogram_intersection(new_helix.T[i], new_helix.T[j])
        mat[j,i] = mat[i,j]
        mat[i,i] = 1
mat[0,0] = 1
plt.imshow(mat)
plt.show()

df = pd.DataFrame(mat)
corrplot(df)#, segid=i+1)
plt.show()


"""
hair_cnt = []
for col in range(hair_arr[0,:].size):
    hair_cnt.append(sum(hair_arr[:,col]))
    
helix_cnt = []
for col in range(helix_arr[0,:].size):
    helix_cnt.append(sum(helix_arr[:,col]))
    
# MINIMUM VALUE > 0
helix_minval = []
for row in range(len(helix_d.array.T)):
    try: helix_minval.append(min(i for i in helix_d.array.T[row,1:] if i > 0))
    except: helix_minval.append(0)

hair_minval = []
for row in range(len(hair_d.array.T)):
    try: hair_minval.append(min(i for i in hair_d.array.T[row,1:] if i > 0))
    except: hair_minval.append(0)

# COMBINED COUNT
helix_ccnt = []
for row in range(len(helix_d.array.T)):
    helix_ccnt.append(sum(helix_d.array.T[row,1:]))

hair_ccnt = []
for row in range(len(hair_d.array.T)):
    hair_ccnt.append(sum(hair_d.array.T[row,1:]))
    

# SECOND MAX
helix_smax = []
for row in range(len(helix_d.array.T)):
    helix_smax.append(sorted(helix_d.array.T[row,1:])[-2])

hair_smax = []
for row in range(len(hair_d.array.T)):
    hair_smax.append(sorted(hair_d.array.T[row,1:])[-2])

# RELATIVE NUMBER OF SAMPLES

helix_ncnt = 0
for row in range(len(helix_d.array.T)):
    helix_ncnt += sum(helix_d.array.T[row,1:])
    
hair_ncnt = 0
for row in range(len(hair_d.array.T)):
    hair_ncnt += sum(hair_d.array.T[row,1:])
    
helix_rel = helix_ccnt/helix_ncnt
hair_rel = hair_ccnt/hair_ncnt

# PLOTTING

plt.plot(helix_d.array[0], helix_minval, label = 'Helix min')
plt.plot(hair_d.array[0], hair_minval, label = 'Hair min')

plt.plot(helix_d.array[0], helix_ccnt, label = 'Helix ccnt')
plt.plot(hair_d.array[0], hair_ccnt, label = 'Hair ccnt')    

plt.plot(helix_d.array[0], helix_smax, label = 'Helix smax')    
plt.plot(hair_d.array[0], hair_smax, label = 'Hair smax')

plt.plot(helix_d.array[0], helix_rel, label = 'Helix rel')    
plt.plot(hair_d.array[0], hair_rel, label = 'Hair rel')

plt.legend()
plt.show()


# PERCENTAGE OVERLAP

helix_max = []
for row in range(len(helix_d.array.T)):
    helix_max.append(max(helix_d.array.T[row,1:]))
hair_max = []
for row in range(len(hair_d.array.T)):
    hair_max.append(sorted(hair_d.array.T[row,1:]))
"""
