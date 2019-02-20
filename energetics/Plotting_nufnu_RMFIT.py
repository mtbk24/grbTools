'''

This program is set up to read the json files written during our PyXSPEC runs. 
ex: 
	filename = '/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/%s/%s/xspec_modelreload_%s_%s_%s_.json'%(burst, detdir, model, model, version, det)
 
The results from each PyXspec run are saved in these files as a dictionary of model components and parameters/settings, 
our model functions in this file are set up to be compatible. 
	For example, the cutoffpl model function uses parameter names 'PhoIndex' and 'HighECut'. 


The Calc_Flux_FromFile(emin, emax, redshift, model, **Pars) function calls the appropriate model functions 
and integrates them over emin to emax. 


We do not need to import any model functions since we have them within this file. 

This program is designed to integrate one model function at a time.  For example, the model 'grbm+bbody+lpow'
	must be passed to Calc_Flux_FromFile three times, once for grbm, once for bbody, and once for lpow. We can compute
	individual flux contributions of each model to the final flux, as a percentage. 




'''




from __future__ import division
import os
import numpy as np
import pandas as pd
from scipy import integrate
from math import pi
from numpy import exp, power, float64, array, inf, logical_and, log, log10, exp
from math import atanh
from scipy.integrate import quad, quadrature
import json
import matplotlib.pyplot as plt
from collections import OrderedDict

from Zoldak.Plotting.publications import journal, journal2
journal2()

data = pd.DataFrame.from_csv(path='/Users/KimiZ/GRBs/analysis/Results/RMFIT_results.csv', header=0, sep=',', index_col=None)








def grbm(energy, **parsdict):
	'''
		RMFIT's BAND function.


	'''
	eng	 	= energy
	a 		= parsdict['alpha']
	b 		= parsdict['beta']
	Ep		= parsdict['epeak']
	N 		= parsdict['norm']
	cond1 = eng < (a-b)*Ep/(2+a)
	cond2 = eng >= (a-b)*Ep/(2+a)
	ans = np.piecewise(eng, [cond1, cond2],\
		[lambda x: N*(power(x/100., a) * exp(-x*(2+a)/Ep) ), \
		 lambda x: N*(power((a-b)*Ep/(100.*(2+a)),a-b)*exp(b-a)*power(x/100.,b))])
	return ans



def sbpl(energy, **parsdict):
	'''
	energy: 	array, energy array for integration. 
	alpha: 		float, low-energy index
	beta: 		float, high-energy index
	enterm: 	float, characteristic energy of function
	norm: 		float, model normalization or amplitude. MUST BE UNLOGGED!!
	entype: 	str, 'epeak', 'E0', or 'ebreak'. MUST ALWAYS SPECIFY NOW. 

	Three different ways to use this function. 
	------------------------------------------
	#1
	Parameter dictionary, if parsdict is used.
	p = {'alpha':-1.25, 
	 		'beta':-2.55, 
	 		'enterm': 300.4, 
	 		'norm': 0.006, 
	 		'entype': 'ebreak'}
	SBPL(energy, **p)

	#2
	Parameter array, if pars is used
	p = [-1.25, -2.55, 300.4, 0.006, 'ebreak']
	SBPL(energy, *p)

	#3
	If parameters are entered individually, in order, with *pars and **parsdict ignored.
	SBPL(energy, -1.25, -2.55, 300.4, 0.006, 'ebreak')

	------------------------------------------
	------------------------------------------

		Conversion between epeak and ebreak:
	epeak = (ebreak) * (10**(0.3 * atanh((a+b+4)/(a-b))))
	ebreak = (epeak) / (10**(0.3 * atanh((a+b+4)/(a-b))))
	# where a and b are alpha and beta

	'''
	a 			= parsdict['alpha']
	b 			= parsdict['beta']
	k 			= parsdict['ebreak']
	N 			= parsdict['norm']
	eng  		= energy
	d			= 0.3 
	# SEVERAL CONSTANTS USED MULTIPLE TIMES. 
	p1   = (b-a)/2.
	p2   = (a+b)/2.
	p3   = ( log10(100.0/k)/d )
	p4   = ( log10(eng/k)/d )
	return N * ((eng/100.0)**p2) * (10**((p1 * d * log((exp(p4) + exp(-p4))/2.)) - (p1 * d * log((exp(p3) + exp(-p3))/2.))))


def cutoffpl(energy, **parsdict):
	'''
			XSPEC's CUTOFFPL function.
			If using function N * ((eng/100.0)**a) * (exp(-eng/t)), you must 
	FLIP SIGN ON PHOTON INDEX AND MULTIPLY THE NORMALIZATION BY 100^(-PHOTON INDEX)
			
	
	'''
	a 			= parsdict['alpha']
	Ep 			= parsdict['epeak']
	N 			= parsdict['norm']
	eng 		= energy
	#return (N*(100**(-a))) * ((eng/100.0)**(-a)) * (exp(-eng/t))
	return N*exp(-eng*(2+a)/Ep )* power(eng/100.0, a)
	#return N * ((eng/100.0)**a) * (exp(-eng/t))



def bbody(energy, **parsdict):
	'''
	XSPEC Blackbody

	------------------------------------------
	------------------------------------------

	This function is: 
		N * (((eng**2)*(8.0525)) / ((kT**4) * (exp(eng/kT)-1)))

	Normalized differently from RMFIT's.
	'''
	kT 			= parsdict['kT']
	N 			= parsdict['norm']
	eng  		= energy
	return N*power(eng,2)*power(exp(eng/float64(kT))-1,-1)
	#return N * (((eng**2)*(8.0525)) / ((kT**4) * (exp(eng/kT)-1)))


def lpow(energy, **parsdict):
	'''
		RMFIT Power-law.


	'''
	a 			= parsdict['alpha']
	N 			= parsdict['norm']
	eng  		= energy
	return N * ((eng/100.0)**a)





def ModelComponentFunc(model, **parameters):
	'''
	
	'''

	listofmodels = ['band', 'sbpl', 'copl']

	if model == 'band':
		function = lambda energy: grbm(energy, **parameters)
	
	elif model == 'sbpl':
		function = lambda energy: sbpl(energy, **parameters)
	
	elif model == 'copl':
		function = lambda energy: cutoffpl(energy, **parameters)
	
	else:
		raise Exception, "Don't recognize models. Options are: %r"%listofmodels
	return function




def get_parameters(burst, model, detector):
	d = data[(data.trigger == burst) & (data.instrument == detector) & (data.function == model)]
	pars = OrderedDict()
	if model == 'band':
		pars['alpha'] = float(d.alpha)
		pars['beta'] = float(d.beta)
		pars['epeak'] = float(d.epeak)
		pars['norm'] = float(d.norm)
	
	elif model == 'sbpl':
		pars['alpha'] = float(d.alpha)
		pars['beta'] = float(d.beta)
		pars['ebreak'] = float(d.ebreak)
		pars['norm'] = float(d.norm)
	
	elif model == 'copl':
		pars['alpha'] = float(d.alpha)
		pars['beta'] = float(d.beta)
		pars['epeak'] = float(d.epeak)
		pars['norm'] = float(d.norm)
	else:
		raise Exception, "Don't recognize models. Options are: %r"%listofmodels
	return pars




def PLT_NuFnu(model, nufnutype, ylims=None, ax=None, **pltArgs): #**parameters):
	'''
	nufnutype= 'individual', 'combined', or 'both'. 

	'''

	if ax is None:
		ax = plt.gca()

	# ENERGY ARRAY, WILL BE DATA ALONG X-AXIS.
	x 	= np.logspace(0.01, 8, 4000)
	# GET PARAMETERS IF DON'T HAVE ANY.
	#if not parameters:
	#parameters = get_parameters(burst='bn100414097', model='band', detector='GBM.fit')
	parameters = get_parameters(burst=burst, model=model, detector='GBM.fit')
	# SPLIT MODEL UP INTO ITS INDIVIDUAL COMPONENTS TO GET EACH CONTRIBUTION. WILL ADD TOGETHER AT END. 
	parts = []
	for mod in model.split('+'):
		parts.append(x * x * ModelComponentFunc(mod, **parameters)(x))

	if ('combined' in nufnutype) or ('both' in nufnutype):
		PLT = ax.plot(x, sum(parts), **pltArgs) # color='blue', ls='-', lw=2, alpha=0.5)
	if ('individual' in nufnutype) or ('both' in nufnutype):
		for i in parts:
			ax.plot(x, i, color='grey', ls='--', lw=1, alpha=0.7)	   
	ax.set_xlim(1.0, 3E7)
	if ylims is not None:
		ax.set_ylim(ylims[0], ylims[1])
	ax.set_xscale('log')
	ax.set_yscale('log')
	return PLT


def PLT_IntegrationRange(redshift, lower=None, upper=None, fillbetween=True, ax=None):
	'''
	lower, upper -- integrand energies. 
	lower = 10, upper = 10,000
	fillbetween

	'''
	if ax is None:
		ax = plt.gca()

	if (lower is None) or (upper is None):
		lower = 10.
		upper = 10000.

	if fillbetween is True:
		PLT = ax.fill_between([10**(-14), lower/(1.+redshift)], 10**(-14), 10**14, 
								color='grey', alpha=0.2)
		ax.fill_between([upper/(1.+redshift), 10**8], 10**(-14), 10**14, 
								color='grey', alpha=0.2)
	else:
		PLT = ax.vlines([10/(1.+redshift), 1E4/(1.+redshift), 1E7/(1.+redshift)], 10**(-14), 10**14, 
					color='grey', linestyle='-.', alpha=0.9)
	return PLT


def PLT_DetectorCoverage(fillbetween=True, ax=None):
	if ax is None:
		ax = plt.gca()

	if fillbetween is True:
		PLT = ax.fill_between([8, 38000], 10**(-14), 10**14, 
								color='pink', alpha=0.1)
		ax.fill_between([2E4, 3E8], 10**(-14), 10**14, 
								color='cyan', alpha=0.1)
	else:
		PLT = ax.hlines(890, 8, 38000, color='pink', lw=6, alpha=0.3)	   # GBM
		ax.hlines(600, 2E4, 3E8, color='goldenrod', lw=6, alpha=0.3)		# LAT			
	return PLT







burst   	= 'bn100414097'
z 			= 1.368
yLims 		= [1E-1, 5E3]
models 		= ['band', 'sbpl', 'copl']  
dets 		= ['GBM.fit', ]*3


# burst   	= 'bn120711115'
# z 			= 1.41
# yLims 		= [1E0, 5E3]
# models 		= ['band', 'sbpl', 'copl']  
# dets 		= ['GBM.fit', ]*3



colordict = {'band':'red',
			 'sbpl':'dodgerblue',
			'copl':'cyan'}
			# 'grbm+lpow':'magenta',
			# 'grbm+bbody':'maroon',
			# 'grbm+bbody+lpow':'brown',
			# 'sbpl+lpow':'purple',
			# 'sbpl+bbody':'navy',
			# 'sbpl+bbody+lpow':'black',
			# 'cutoffpl+lpow':'lime',
			# 'cutoffpl+bbody':'green',
			# 'cutoffpl+bbody+lpow':'forestgreen'}



plt.clf()
PLT_DetectorCoverage(fillbetween=True, ax=None)
PLT_IntegrationRange(redshift=z, lower=10., upper=1E7, fillbetween=False, ax=None)

for model,det in zip(models,dets):
	PLT_NuFnu(model, nufnutype='combined', ylims=yLims, ax=None, 
				**dict(color=colordict[model], ls='--', lw=1.5, alpha=0.5, label=model))
plt.legend(loc=8, numpoints=1, labelspacing=0.1, handletextpad=0.2)
plt.xlabel('Energy (keV)')
plt.ylabel(r'$\nu F_{\nu}$')
plt.tight_layout(pad=1)
outfilename = '/Users/KimiZ/GRBs2/python_modeling/%s_nufnu_severalmodels.pdf'%burst
plt.savefig(outfilename)
os.system('open %s'%outfilename)

#plt.show()













