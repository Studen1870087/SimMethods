#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 15:13:40 2019

@author: jake
"""

import numpy as np
import matplotlib.pyplot as plt

"LOAD DATA"

prot = 'hairpin'
file = '/force-vs-distance'
ext = '.dat'

data = np.loadtxt( prot+file+ext )

"PLOT"
plt.plot(data[:,0], data[:,1])