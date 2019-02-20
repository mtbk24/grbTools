
from __future__ import division
import os, glob
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import spectralTools
from spectralTools.models import BAND, SBPL, COPL, LPOW, BBODY

from GRBAnalysis.analysis import Data
from Zoldak.Cosmology.luminositydistance import LumDist
from Zoldak.Flux.calculateFluxFlueEiso import Calc_Eiso, Calc_Flux


'''

BAND(energy, alpha, beta, enterm, norm, entype='E0')
SBPL(energy, alpha, beta, enterm, norm, entype='ebreak')
COPL(energy, index, enterm, norm, entype='E0')

LPOW(energy, index, norm)
BBODY(energy, kT, norm)


'''


def get_nPars(modelname):
	if modelname in ('grbm', 'sbpl'):
		return 4
	elif modelname in ('lpow', 'bbody'):
		return 2
	elif modelname in ('cutoffpl'):
		return 3
	else:
		raise Exception, "Don't recognize model"


def get_modeldict(modelname, pars):
	if modelname == 'grbm':
		return {'grbm': {
                        'alpha':    pars[0],
                        'beta': 	pars[1], 
                        'enterm':   pars[2],
                        'norm':     pars[3],
                        'entype':   'E0'
                        }}
	elif modelname == 'sbpl':
		return {'sbpl': {
                        'alpha':    pars[0],
                        'beta': 	pars[1], 
                        'enterm':   pars[2],
                        'norm':     pars[3],
                        'entype':   'ebreak'
                        }}
	elif modelname == 'cutoffpl':
		return {'cutoffpl': {
                        'index':    pars[0],
                        'enterm':   pars[1],
                        'norm':     pars[2],
                        'entype':   'E0'
                        }}
	elif modelname == 'lpow':
		return {'lpow': {
                        'index':    pars[0],
                        'norm':     pars[1],
                        }}
	elif modelname == 'bbody':
		return {'bbody': {
                        'kT':    	pars[0],
                        'norm':     pars[1],
                        }}
	else:
		raise Exception, "Don't recognize model"


from __future__ import division
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
from scipy import integrate
from math import pi
from Zoldak.Cosmology.luminositydistance import LumDist

import spectralTools

from spectralTools.models import BAND, SBPL, COPL, BBODY, LPOW



def Calc_Flux(model, Pars, emin=10., emax=10000., redshift=None):
    '''
    Calc_Flux(model, Pars, emin=10., emax=10000., redshift=None)
    
    model:  str, model name.
            ex:  'grbm', 'grbm+bbody', 'grbm+bbody+lpow'
            all individual model options:
                grbm, sbpl, cutoffpl, lpow, and bbody
    
    Pars:   dict, model parameters.  SEE BELOW.
    
    emin:   float, lower energy integration limit. Do not k-correct yourself.
    
    emax:   float, upper energy integration limit. Do not k-correct yourself.
    
    redshift: float, redshift of GRB for k-correcting. If not k-correcting, do not enter anything.
    
    kcorrect:   [True|False]  
                True will use k-correction: emin/(1.+redshift) to emax/(1.+redshift)
    
    
    
    # REGARDLESS OF WHETHER THERE IS ONE MODEL OR MORE THAN ONE, PARAMETER DICTIONARIES SHOULD
    # ALWAYS BE SET UP THIS WAY.  THIS IS TO ENSURE PARAMETERS BETWEEN TWO DIFFERENT MODELS DO NOT GET MIXED UP.
    
    pars = {'grbm':     {'alpha':   -1.0,
                         'beta':    -2.5,
                         'enterm':  300.0,
                         'norm':    0.01,
                         'entype':  'epeak'},
            'sbpl':     {'alpha':   -1.0,
                         'beta':    -2.5,
                         'enterm':  300.0,
                         'norm':    0.01,
                         'entype':  'epeak'},
            'cutoffpl': {'alpha':   -1.0,
                         'enterm':  300.0,
                         'norm':    0.01,
                         'entype':  'epeak'},
            'lpow':     {'index':   -1.0,
                         'norm':    0.001},
            'bbody':    {'kT':      10.0,
                         'norm':    0.001} } 
                     
    
    '''
    keVtoerg    = 1.60217657E-9
    if redshift is not None:
        emin_K      = emin/(1.+redshift)
        emax_K      = emax/(1.+redshift)
    else:
        print('\n *** WARNING:  NOT USING K-CORRECTED ENERGIES. *** \n')
        emin_K      = emin
        emax_K      = emax

    if 'grbm' == model:
        function = lambda energy: energy * BAND(energy,  
                                                P['grbm']['alpha'], 
                                                P['grbm']['beta'], 
                                                P['grbm']['enterm'], 
                                                P['grbm']['norm'], 
                                                entype=P['grbm']['entype'])
    
    if 'sbpl' == model:
        function = lambda energy: energy * SBPL(energy,  
                                    P['sbpl']['alpha'], 
                                    P['sbpl']['beta'], 
                                    P['sbpl']['enterm'], 
                                    P['sbpl']['norm'], 
                                    entype=P['sbpl']['entype'])
    
    if 'cutoffpl' == model:
        function = lambda energy: energy * COPL(energy,  
                                    P['cutoffpl']['alpha'], 
                                    P['cutoffpl']['enterm'], 
                                    P['cutoffpl']['norm'], 
                                    entype=P['cutoffpl']['entype'])
                                                
                                                
    
    if 'lpow' == model:
        function = lambda energy: energy * LPOW(energy,  
                                    P['lpow']['alpha'], 
                                    P['lpow']['norm'])
                                                
    
    
    if 'bbody' == model:
        function = lambda energy: energy * BBODY(energy,  
                        P['bbody']['kT'], 
                        P['bbody']['norm'])

    Flux = integrate.quad(function, emin_K, emax_K, limit=100)[0] * keVtoerg
    nrgFlux.append(Flux)
    nrgFlux      = sum(nrgFlux)
    return nrgFlux






burst 		= 'bn080916009'
detector   	= 'L'
#model 		= 'grbm+bbody+lpow'
model = 'grbm'
version 	= '-01-'




#def run_program():
params_dict = {
				'burst':        burst,
				'det':          detector,
				'model':        model,
				'pyx_version':  version}

test    		= Data(**params_dict)   
redshift        = test._redshift
duration	    = test._t90stop - test._t90start
parameters      = test.pyx_bestfit()['parameters']   # parameters


modparts = model.split('+')
nParts = len(modparts)


AllFluxes = []
AllEisos = []
n = 0
for mod in modparts:
	npars 		= get_nPars(mod)
	print(n, n+npars)
	p 		= parameters[n : n+npars] # SUBSET OF PARS BASED ON MOD
	pardict 	= get_modeldict(mod, p)
	Flux 		= Calc_Flux(mod, pardict, 10., 10000000., redshift) 
	Eiso        = Calc_Eiso(mod, pardict, 10., 10000000., duration, redshift) 
	AllFluxes.append(Flux)
	AllEisos.append(Eiso)
	n = n+npars




# for i in range(0, nParts):
# 	print(i)







# pars1 = {'grbm': {
#                         'alpha':    p[0],
#                         'beta': 	p[1], 
#                         'enterm':   p[2],
#                         'norm':     p[3],
#                         'entype':   'E0'
#                         }}
# pars2 = {'lpow': {
#                         'alpha':    p[4], 
#                         'norm':     p[5],
#                         }}


# # Eiso1       = Calc_Eiso('grbm', pars1, 10., 10000000., duration, redshift) 
# # Flux1 		= Calc_Flux('grbm', pars1, 10., 10000000., redshift) 
# # Eiso2       = Calc_Eiso('lpow', pars2, 10., 10000000., duration, redshift) 
# # Flux2 		= Calc_Flux('lpow', pars2, 10., 10000000., redshift) 

# Eiso1       = Calc_Eiso('grbm', pars1, 10., 10000000., duration, redshift) 
# Flux1 		= Calc_Flux('grbm', pars1, 10., 10000000., redshift) 
# Eiso2       = Calc_Eiso('lpow', pars2, 10., 10000000., duration, redshift) 
# Flux2 		= Calc_Flux('lpow', pars2, 10., 10000000., redshift) 

# print('%9.3E %9.3E %9.3E  6.2%f %6.2f'%(Flux1, Flux2, Flux1+Flux2, Flux1/(Flux1+Flux2), Flux2/(Flux1+Flux2)))
# print('%9.3E %9.3E %9.3E  6.2%f %6.2f'%(Eiso1, Eiso2, Eiso1+Eiso2, Eiso1/(Eiso1+Eiso2), Eiso2/(Eiso1+Eiso2)))









# params_dict = {
# 				'burst':        'bn80916009',
# 				'det':          'L',
# 				'model':        'grbm',
# 				'pyx_version': '-01-',}

# test    = Data(**params_dict)   
# z       = test._redshift
# dur     = test._t90stop - test._t90start
# p       = test.pyx_bestfit()['parameters'] # parameters


# pars1 = {'grbm': {
#                         'alpha':    p[0],
#                         'beta': 	p[1], 
#                         'enterm':   p[2],
#                         'norm':     p[3],
#                         'entype':   'E0'
#                         }}

# Eiso1       = Calc_Eiso('grbm', pars1, 10., 10000000., dur, z) 




