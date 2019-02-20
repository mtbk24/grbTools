


from __future__ import division
import os, glob
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import spectralTools
from spectralTools import models as mod

from GRBAnalysis.analysis import Data

burst       = 'bn100728095'
model       = 'grbm+bbody'
version     = '-01-'


# burst       = 'bn091208410'
# model       = 'cutoffpl+bbody'
# version     = '-01-'

det         = 'L'
detdir      = 'GBMwLAT' if 'L' in det else 'GBM'

params_dict = {
'burst':        burst,
'det':          det,
'model':        model,
'pyx_version': version,}

test            = Data(**params_dict)   
z               = test._redshift
t90start        = test._t90start
t90stop         = test._t90stop
duration        = t90stop - t90start
parameters      = test.pyx_bestfit()['parameters'] # parameters
p               = parameters
energy_limits   = [10.0, 10000.0]
eng_lim_kcor    = [i/(1.0+z) for i in energy_limits]


x = np.logspace(0.01, 8, 4000)


y1 = mod.COPL(energy = x, 
             index  = p[0], 
             enterm = p[1], 
             norm   = p[2], 
             entype = 'E0')
#mod.BBODY(energy, kT, N)
y2 = mod.BBODY(energy = x, 
             kT  = p[3],  
             norm   = p[4]) 

Y = y1+y2


ylims = [1E0, 10**3.2]


plt.figure(figsize=(7,4))
plt.grid(which='major', axis='both', alpha=0.2)
plt.fill_between([.001, 10./(1.+z)], 10**(-14), 10**4, 
                  color='grey', alpha=0.2)
plt.fill_between([10000./(1.+z), 10**8], 10**(-14), 
                  10**4, color='grey', alpha=0.2)

plt.plot(x, x*x*y1, color='grey', ls='--', lw=1, alpha=0.7)
plt.plot(x, x*x*y2, color='green', ls='--',lw=1, alpha=0.7)
plt.plot(x, x*x*Y, color='blue', ls='-', lw=1, alpha=0.7)
plt.hlines(890, 8, 38000, color='pink', lw=6, alpha=0.3)        # GBM
plt.hlines(600, 2E4, 3E8, color='goldenrod', lw=6, alpha=0.3)        # LAT            
plt.loglog()

plt.xlim(1.0, 10**8)
plt.ylim(ylims[0], ylims[1])
#plt.ylim(1E-2, 10**3.2)
plt.xlabel('Energy (keV)', fontsize=12)
plt.ylabel(r'$\nu F_{\nu}$', fontsize=16)
plt.tight_layout(pad=0.1, w_pad=0.1, h_pad=0.1) 
plt.show()












