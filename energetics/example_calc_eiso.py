'''
This program is an example of how to calculate Eiso energies using the Data saved from PYXSPEC fits. 


# params_dict = {
# 'burst':        burst,
# 'det':          det,
# 'model':        model,
# 'pyx_version': version,
# 'bxa_version': '-01-',
# 'bxa_Aversion':'-01-'
# }


'''


from __future__ import division
from collections import OrderedDict
import os, glob
import numpy

from GRBAnalysis.analysis import Data
from Zoldak.Cosmology.luminositydistance import LumDist
from Zoldak.Flux.calculateFluxFlueEiso import Calc_Eiso, Calc_Flux
from Zoldak.Math.tools import percentage_difference


params_dict = {
				'burst':        'bn090902462',
				'det':          'L',
				'model':        'cutoffpl+bbody+lpow',
				'pyx_version': '-01-',}



test    = Data(**params_dict)   
z       = test._redshift
dur     = test._t90stop - test._t90start
p       = test.pyx_bestfit()['parameters'] # parameters


pars1 = {'cutoffpl': {
                        'alpha':    p[0], 
                        'enterm':   p[1],
                        'norm':     p[2],
                        'entype':   'E0'
                        }}

pars2 = {'bbody': {
                        'kT':     p[3], 
                        'norm':   p[4],
                        }}
pars3 = {'lpow': {
                        'alpha':    p[5], 
                        'norm':     p[6],
                        }}

Eiso1       = Calc_Eiso('cutoffpl', pars1, 10., 10000., dur, z) + \
			 Calc_Eiso('bbody', pars2, 10., 10000., dur, z) + \
			 Calc_Eiso('lpow', pars3, 10., 10000., dur, z)


Eiso2       = Calc_Eiso('cutoffpl', pars1, 10., 10000000., dur, z) + \
			 Calc_Eiso('bbody', pars2, 10., 10000000., dur, z) + \
			 Calc_Eiso('lpow', pars3, 10., 10000000., dur, z)

percentage_difference(Eiso1, Eiso2)










params_dict = {
				'burst':        'bn090902462',
				'det':          'L',
				'model':        'cutoffpl+lpow',
				'pyx_version': '-01-',}

test    = Data(**params_dict)   
z       = test._redshift
dur     = test._t90stop - test._t90start
p       = test.pyx_bestfit()['parameters'] # parameters

pars1 = {'cutoffpl': {
                        'alpha':    p[0], 
                        'enterm':   p[1],
                        'norm':     p[2],
                        'entype':   'E0'
                        }}

pars2 = {'lpow': {
                        'alpha':    p[3], 
                        'norm':     p[4],
                        }}

Eiso1       = Calc_Eiso('cutoffpl', pars1, 10., 10000., dur, z) + \
			 Calc_Eiso('lpow', pars2, 10., 10000., dur, z)


Eiso2       = Calc_Eiso('cutoffpl', pars1, 10., 10000000., dur, z) + \
			 Calc_Eiso('lpow', pars2, 10., 10000000., dur, z)

percentage_difference(Eiso1, Eiso2)














params_dict = {
				'burst':        'bn090902462',
				'det':          'L',
				'model':        'grbm',
				'pyx_version': '-01-',}

test    = Data(**params_dict)   
z       = test._redshift
dur     = test._t90stop - test._t90start
p       = test.pyx_bestfit()['parameters'] # parameters


pars1 = {'grbm': {
                        'alpha':    p[0],
                        'beta': 	p[1], 
                        'enterm':   p[2],
                        'norm':     p[3],
                        'entype':   'E0'
                        }}

Eiso1       = Calc_Eiso('grbm', pars1, 10., 10000., dur, z)

Eiso2       = Calc_Eiso('grbm', pars1, 10., 10000000., dur, z) 

percentage_difference(Eiso1, Eiso2)















LATbursts = ['bn080916009', 'bn090323002', 'bn090328401','bn090510016',
			'bn090902462','bn090926181','bn091003191','bn091208410',
			'bn100728095','bn110731465','bn130518580','bn131108862','bn131231198']



for burst in LATbursts:
	params_dict = {
					'burst':        burst,
					'det':          'L',
					'model':        'grbm',
					'pyx_version': '-01-',}

	test    	= Data(**params_dict)   
	z       	= test._redshift
	duration   	= test._t90stop - test._t90start
	p       	= test.pyx_bestfit()['parameters'] # parameters


	pars1 = {'grbm': {
	                        'alpha':    p[0],
	                        'beta': 	p[1], 
	                        'enterm':   p[2],
	                        'norm':     p[3],
	                        'entype':   'E0'
	                        }}

	Flux 		= Calc_Flux(model='grbm', Pars=pars1, emin=10., emax=10000., redshift=z)
	Fluence 	= Flux * duration
	Eiso        = Calc_Eiso('grbm', pars1, 10., 10000., duration, z)
	#print('%s:  %9.3E  %9.3E  %9.3E'%(burst, Flux, Fluence, Eiso))
	#print('%s,%.3E,%.3E,%.3E'%(burst, Flux, Fluence, Eiso))
	print('%.2f'%duration)







