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
from scipy import integrate
from math import pi
from numpy import exp, power, float64, array, inf, logical_and, log, log10, exp
from math import atanh
from scipy.integrate import quad, quadrature
import json
import matplotlib.pyplot as plt

from Zoldak.Plotting.publications import journal, journal2
journal2()


def grbm(energy, **parsdict):
	'''
		XSPEC's BAND function.


	'''

	eng     = energy
	a 		= parsdict['alpha'][0]
	b 		= parsdict['beta'][0]
	t		= parsdict['tem'][0]
	N 		= parsdict['norm'][0]

	cond1 = eng < (a-b) * t
	cond2 = eng >= (a-b) * t
	
	return np.piecewise(eng, [cond1, cond2],\
			[lambda x: N * (power(x/100., a) * exp(-x/t)), \
			 lambda x: N * (((((a-b)*t)/100.0)**(a-b)) * (exp(b-a))*((x/100.0)**b))])
#
#	if eng < ( (a-b) * t ):
#		return  N * (((eng/100.0)**a) * (exp(-eng/t)))
#	else:
#		return  N * (((((a-b)*t)/100.0)**(a-b)) * (exp(b-a))*((eng/100.0)**b))




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
	a 			= parsdict['alpha'][0]
	b 			= parsdict['beta'][0]
	k 			= parsdict['ebreak'][0]
	N 			= parsdict['norm'][0]
	eng  		= energy
	d			= 0.3 
	# SEVERAL CONSTANTS USED MULTIPLE TIMES. 
	p1   = (b-a)/2.
	p2   = (a+b)/2.
	p3   = ( log10(100.0/k)/d )
	p4   = ( log10(eng/k)/d )
	return N * ((eng/100.0)**p2) * (10**((p1 * d * log((exp(p4) + exp(-p4))/2.)) - (p1 * d * log((exp(p3) + exp(-p3))/2.))))



def bbody(energy, **parsdict):
	'''
	XSPEC Blackbody

	------------------------------------------
	------------------------------------------

	This function is: 
		N * (((eng**2)*(8.0525)) / ((kT**4) * (exp(eng/kT)-1)))

	Normalized differently from RMFIT's.
	'''
	kT 			= parsdict['kT'][0]
	N 			= parsdict['norm'][0]
	eng  		= energy
	return N * (((eng**2)*(8.0525)) / ((kT**4) * (exp(eng/kT)-1)))

	#if eng <= (709.666 * kT): # to avoid exp overflow error
	#	return N * (((eng**2)*(8.0525)) / ((kT**4) * (exp(eng/kT)-1)))
	#else:
	#	return 0



def lpow(energy, **parsdict):
	'''
		XSPEC Power-law, should be same as the RMFIT's. Might need a sign 
	change on alpha. 


	'''
	a 			= parsdict['plIndex'][0]
	N 			= parsdict['norm'][0]
	eng  		= energy
	return N * ((eng/100.0)**a)



def cutoffpl(energy, **parsdict):
	'''
			XSPEC's CUTOFFPL function.
			If using function N * ((eng/100.0)**a) * (exp(-eng/t)), you must 
	FLIP SIGN ON PHOTON INDEX AND MULTIPLY THE NORMALIZATION BY 100^(-PHOTON INDEX)
			
	
	'''
	a 			= parsdict['PhoIndex'][0]
	t 			= parsdict['HighECut'][0]
	N 			= parsdict['norm'][0]
	eng 		= energy
	return (N*(100**(-a))) * ((eng/100.0)**(-a)) * (exp(-eng/t))
	#return N * ((eng/100.0)**a) * (exp(-eng/t))





def Calc_Flux_FromFile(emin, emax, redshift, model, **Pars):
    '''
    Calc_Flux(model, Pars, emin, emax, redshift=None)
    
    model:  str, model name.
            examples:  
            'grbm', 'sbpl', 'cutoffpl', 'bbody', 'lpow'
            if two-component models are preferred, run this calculation for 
            each individual model and sum the individual components. 
            i.e.:
            
            model = 'grbm+lpow'
            fluxes = []
            for mod in model:
                flux = Calc_Flux(model=mod, Pars=p, emin=10., emax=10000., redshift=1.23)
                fluxes.append(flux)
            sum(fluxes)

    Pars:   dict, model parameters.  **SEE BELOW FOR ALL EXAMPLES.**
    emin:   float, lower energy integration limit. Do not k-correct yourself.
    emax:   float, upper energy integration limit. Do not k-correct yourself.
    redshift: float, redshift of GRB for k-correcting. If not k-correcting, do not enter anything.
    
    
    # PARAMETER DICTIONARIES SHOULD ALWAYS BE SET UP THIS WAY.

	BAND(energy, **p)
    
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
            'lpow':     {'alpha':   -1.0,
                         'norm':    0.001},
            'bbody':    {'kT':      10.0,
                         'norm':    0.001} } 

	
	BAND(energy, **pars['grbm'])
	SBPL(energy, **pars['sbpl'])
	COPL(energy, **pars['cutoffpl])
	LPOW(energy, **pars['lpow'])
	BBODY(energy, **pars['bbody'])
                     
    
    '''

    listofmodels = ['grbm', 'sbpl', 'cutoffpl', 'bbody', 'lpow']

    keVtoerg     = 1.60217657E-9
    if redshift is not None:
        emin      = emin/(1.+redshift)
        emax      = emax/(1.+redshift)
    else:
        print('\n *** WARNING:  NOT USING K-CORRECTED ENERGIES. *** \n')
        pass
    if model == 'grbm':
        function = lambda energy: energy * grbm(energy, **Pars['grbm'])
    
    elif model == 'sbpl':
        function = lambda energy: energy * sbpl(energy, **Pars['sbpl'])
    
    elif model == 'cutoffpl':
        function = lambda energy: energy * cutoffpl(energy, **Pars['cutoffpl'])
                                               
    elif model == 'lpow':
        function = lambda energy: energy * lpow(energy, **Pars['lpow'])
    
    elif model == 'bbody':
        function = lambda energy: energy * bbody(energy, **Pars['bbody'])
    
    else:
        raise Exception, "Don't recognize models. Options are: %r"%listofmodels

    nrgFlux = integrate.quad(function, emin, emax, limit=100)[0] * keVtoerg
    return nrgFlux

   

def ModelComponentFunc(model, **parameters):
    '''
    ModelComponentFunc(model, **parameters)
    
    model:  str, model name.
            examples:  
            'grbm', 'sbpl', 'cutoffpl', 'bbody', 'lpow'
            if two-component models are preferred, run this calculation for 
            each individual model and sum the individual components. 
            i.e.:
            
            model = 'grbm+lpow'
            fluxes = []
            for mod in model:
                flux = Calc_Flux(model=mod, parameters=p, emin=10., emax=10000., redshift=1.23)
                fluxes.append(flux)
            sum(fluxes)

    Pars:   dict, model parameters.  **SEE BELOW FOR ALL EXAMPLES.**

    
    # PARAMETER DICTIONARIES SHOULD ALWAYS BE SET UP THIS WAY.

	BAND(energy, **p)
    
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
            'lpow':     {'alpha':   -1.0,
                         'norm':    0.001},
            'bbody':    {'kT':      10.0,
                         'norm':    0.001} } 

	
	BAND(energy, **pars['grbm'])
	SBPL(energy, **pars['sbpl'])
	COPL(energy, **pars['cutoffpl])
	LPOW(energy, **pars['lpow'])
	BBODY(energy, **pars['bbody'])
                     
    
    '''

    listofmodels = ['grbm', 'sbpl', 'cutoffpl', 'bbody', 'lpow']

    if model == 'grbm':
        function = lambda energy: grbm(energy, **parameters['grbm'])
    
    elif model == 'sbpl':
        function = lambda energy: sbpl(energy, **parameters['sbpl'])
    
    elif model == 'cutoffpl':
        function = lambda energy: cutoffpl(energy, **parameters['cutoffpl'])
                                               
    elif model == 'lpow':
        function = lambda energy: lpow(energy, **parameters['lpow'])
    
    elif model == 'bbody':
        function = lambda energy: bbody(energy, **parameters['bbody'])
    
    else:
        raise Exception, "Don't recognize models. Options are: %r"%listofmodels
    return function



def get_parameters(burst, model, version, detector):
	detdir 		= ('GBM' if 'G' in det else 'GBMwLAT')
	filename 	= '/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/%s/%s/xspec_modelreload_%s_%s_%s_.json'%(burst, detdir, model, model, version, detector)
	pars 		= json.load(open(filename))
	return pars




# def PLT_NuFnu(model, nufnutype, redshift=None, ylims=None, detectorcoverage=None, ax=None, **parameters):
# 	'''
# 	nufnutype= 'individual', 'combined', or 'both'. 

# 	'''

# 	if ax is None:
# 		ax = plt.gca()
# 	if ylims is None:
# 		ylims 	= [None, None]

# 	# ENERGY ARRAY, WILL BE DATA ALONG X-AXIS.
# 	x 	= np.logspace(0.01, 8, 4000)
# 	# GET PARAMETERS IF DON'T HAVE ANY.
# 	if not parameters:
# 		parameters = get_parameters(burst=burst, model=model, version=version, detector=det)
	
# 	# SPLIT MODEL UP INTO ITS INDIVIDUAL COMPONENTS TO GET EACH CONTRIBUTION. WILL ADD TOGETHER AT END. 
# 	parts = []
# 	n=0
# 	for mod in model.split('+'):
# 		parts.append(x * x * ModelComponentFunc(mod, **parameters)(x))
# 		n=n+1

# 	if ('combined' in nufnutype) or ('both' in nufnutype):
# 		PLT = ax.plot(x, sum(modelpts), color='blue', ls='-', lw=2, alpha=0.5)
# 	if ('individual' in nufnutype) or ('both' in nufnutype):
# 		for i in modelpts:
# 			ax.plot(x, i, color='grey', ls='--', lw=1, alpha=0.7)

# 	if redshift:
# 		ax.fill_between([.001, 10./(1.+redshift)], 10**(-14), 10**14, 
# 						color='grey', alpha=0.2)
# 		ax.fill_between([10000./(1.+redshift), 10**8], 10**(-14), 10**14, 
# 						color='grey', alpha=0.2)
	
# 	if detectorcoverage:
# 		ax.hlines(890, 8, 38000, color='pink', lw=6, alpha=0.3)        		# GBM
# 		ax.hlines(600, 2E4, 3E8, color='goldenrod', lw=6, alpha=0.3)        # LAT            
# 	ax.set_xlim(1.0, 10**8)
# 	if ylims:
# 		ax.set_ylim(ylims[0], ylims[1])
# 	ax.set_xscale('log')
# 	ax.set_yscale('log')
# 	ax.set_xlabel('Energy (keV)', fontsize=12)
# 	ax.set_ylabel(r'$\nu F_{\nu}$', fontsize=16)
# 	return PLT




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
	parameters = get_parameters(burst=burst, model=model, version=version, detector=det)
	
	# SPLIT MODEL UP INTO ITS INDIVIDUAL COMPONENTS TO GET EACH CONTRIBUTION. WILL ADD TOGETHER AT END. 
	parts = []
	for mod in model.split('+'):
		parts.append(x * x * ModelComponentFunc(mod, **parameters)(x))

	if ('combined' in nufnutype) or ('both' in nufnutype):
		PLT = ax.plot(x, sum(parts), **pltArgs) # color='blue', ls='-', lw=2, alpha=0.5)
	if ('individual' in nufnutype) or ('both' in nufnutype):
		for i in parts:
			ax.plot(x, i, color='grey', ls='--', lw=1, alpha=0.7)       
	#ax.set_xlim(1.0, 10**8)
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
		#PLT = ax.vlines([lower/(1.+redshift), upper/(1.+redshift)], 10**(-14), 10**14, 
		#		color='grey', linestyle='-.', alpha=0.9)
		PLT = ax.vlines([10/(1.+redshift), 1E4/(1.+redshift), 1E7/(1.+redshift)], 10**(-14), 10**14, 
					color='grey', linestyle='-.', alpha=0.9)
		# ax.text(10/(1.+redshift), ylims[1], r'$\frac{10\text{keV}{(1+z)}$', verticalalignment='top')
		# ax.text(1E4/(1.+redshift), ylims[1], r'$\frac{10 MeV}{(1+z)}$', verticalalignment='top')
		# ax.text(1E7/(1.+redshift), ylims[1], r'$\frac{10 GeV}{(1+z)}$', verticalalignment='top')
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
		PLT = ax.hlines(890, 8, 38000, color='pink', lw=6, alpha=0.3)       # GBM
		ax.hlines(600, 2E4, 3E8, color='goldenrod', lw=6, alpha=0.3)        # LAT            
	return PLT




'''
	PLOT MULTIPLE NUFNU SPECTRUMS ON THE SAME PLOT. CAN INTERCHANGE GBM+LAT ('L') AND GBM ('G') RESULTS.
	GBM PLOTS WILL BECOME DASHED LINES ABOVE THE GBM DETECTOR'S COVERAGE, SHOWING THE 
	PATH OF THE MODEL, BUT THAT THE DATA TRUNCATES.  


'''




# burst 		= 'bn090902462'
# z 			= 1.82
# models 		= ['cutoffpl+bbody+lpow','grbm', 'sbpl+bbody','grbm', 'cutoffpl+bbody+lpow']
# versions 	= ['-01-', '-01-', '-01-', '-01-', '-01-']
# dets 		= ['L', 'L', 'G', 'G', 'G']


# burst 		= 'bn090510016'
# z 			= 0.903
# models 		= ['sbpl+lpow','grbm+lpow', 'sbpl+bbody+lpow','grbm+bbody+lpow', 
# 				'cutoffpl+lpow','cutoffpl+bbody+lpow', 'grbm']
# versions 	= ['-01-',]*len(models) #'-01-', '-01-', '-01-', '-01-']
# dets 		= ['L', ]*len(models) #'L', 'G', 'G', 'G']


# burst 		= 'bn090510016'
# z 			= 0.903
# models 		= ['sbpl+lpow','sbpl+bbody+lpow', 'grbm', 'cutoffpl', 'sbpl', 'grbm']#'sbpl+bbody+lpow','grbm+bbody+lpow', 
# 				#'cutoffpl+lpow','cutoffpl+bbody+lpow', 'grbm']
# versions 	= ['-01-',]*len(models) #'-01-', '-01-', '-01-', '-01-']
# dets 		= ['L', ]*3 + ['G', ]*3  #'L', 'G', 'G', 'G']


# burst 		= 'bn080916009'
# z 			= 4.35
# yLims 		= [2E0, .5E3]
# models 		= ['grbm', 'grbm+lpow', 'sbpl', 'grbm+bbody', 'sbpl+bbody', 
# 			   'grbm', 'cutoffpl+bbody', 'sbpl+bbody']
# versions 	= ['-01-','-01-', '-01-', '-02-', '-02-',
# 			   '-01-', '-02-', '-02-']
# dets 		= ['L', ]*5 + ['G', ]*3


# burst 		= 'bn090902462'
# z 			= 1.82
# yLims 		= [1E1, 0.7E4]
# models 		= ['grbm', 'cutoffpl+lpow', 'cutoffpl+bbody+lpow', 
# 			   'grbm', 'sbpl+bbody', 'cutoffpl+bbody+lpow']
# versions 	= ['-01-','-01-', '-01-', 
# 			   '-01-', '-01-', '-01-']
# dets 		= ['L', ]*3 + ['G', ]*3


# burst 		= 'bn130518580'
# z 			= 2.49
# yLims 		= [2E0, 8E2]
# models 		= ['grbm', 'sbpl','sbpl+bbody', 
# 			   'grbm', 'sbpl', 'sbpl+bbody', 'cutoffpl+bbody']
# versions 	= ['-01-', '-01-', '-01-',
# 			   '-01-', '-01-', '-01-', '-01-']
# dets 		= ['L', ]*3 + ['G', ]*4

# models 		= ['grbm', 'sbpl+bbody', 
# 			   'grbm', 'cutoffpl+bbody', 'sbpl+bbody']
# versions 	= ['-01-', '-01-', 
# 			   '-01-', '-01-', '-01-']
# dets 		= ['L', ]*2 + ['G', ]*3



# burst 		= 'bn110731465'
# z 			= 2.83
# yLims 		= [6E0, 2E3]
# models 		= ['grbm', 'grbm+bbody','sbpl+lpow', 'cutoffpl+bbody+lpow',
# 			   'grbm','grbm+bbody']
# versions 	= ['-01-', '-01-', '-01-','-01-',
# 			   '-01-', '-01-']
# dets 		= ['L', ]*4 + ['G', ]*2



# burst 		= 'bn100728095'
# z 			= 1.57
# yLims 		= [1E0, 3E2]
# models 		= ['grbm', 'sbpl','cutoffpl+bbody',
# 			   'grbm', 'sbpl','cutoffpl+bbody']
# versions 	= ['-01-', '-01-', '-01-',
# 			   '-01-', '-01-','-01-']
# dets 		= ['L', ]*3 + ['G', ]*3


# burst   	= 'bn131231198'
# z 			= 0.642
# yLims 		= [1E1, 2E3]
# models 		= ['grbm', 'sbpl+bbody','grbm+bbody',
# 			   'grbm', 'sbpl+bbody','grbm+bbody']
# versions 	= ['-01-', '-01-', '-01-',
# 			   '-01-', '-01-','-01-']
# dets 		= ['L', ]*3 + ['G', ]*3


# burst   = 'bn090926181'
# z 			= 2.11
# yLims 		= [1E1, 5E3]
# models 		= ['grbm', 'sbpl','grbm+lpow', 'grbm+bbody+lpow',
# 			   'grbm', 'sbpl+bbody', 'grbm+bbody+lpow']
# versions 	= ['-01-', '-01-', '-01-','-01-',
# 			   '-01-', '-01-','-01-']
# dets 		= ['L', ]*4 + ['G', ]*3


# burst   	= 'bn131108862'
# z 			= 2.4
# yLims 		= [1E0, 5E2]
# models 		= ['grbm', 'sbpl', 'sbpl+bbody', 'sbpl+lpow',
# 			   'grbm', 'sbpl', 'cutoffpl']
# versions 	= ['-01-', '-01-', '-01-', '-01-',
# 			   '-01-', '-01-','-01-']
# dets 		= ['L', ]*4 + ['G', ]*3


# burst   	= 'bn090328401'
# z 			= 0.74
# yLims 		= [1E0, 5E2]
# models 		= ['grbm', 'sbpl', 'cutoffpl+lpow', 
# 			   'grbm', 'sbpl', 'cutoffpl']
# versions 	= ['-01-', '-01-', '-01-',  
# 			   '-01-', '-01-','-01-']
# dets 		= ['L', ]*3 + ['G', ]*3


burst   	= 'bn091208410'
z 			= 1.06
yLims 		= [1E0, 5E2]
models 		= ['grbm', 'sbpl', 'grbm+bbody',  
			   'grbm', 'sbpl', 'cutoffpl']
versions 	= ['-01-', '-01-', '-01-',
			   '-01-', '-01-','-01-']
dets 		= ['L', ]*3 + ['G', ]*3




colordict = {'grbm':'red',
			 'sbpl':'dodgerblue',
			'cutoffpl':'cyan',
			'grbm+lpow':'magenta',
			'grbm+bbody':'maroon',
			'grbm+bbody+lpow':'brown',
			'sbpl+lpow':'purple',
			'sbpl+bbody':'navy',
			'sbpl+bbody+lpow':'black',
			'cutoffpl+lpow':'lime',
			'cutoffpl+bbody':'green',
			'cutoffpl+bbody+lpow':'forestgreen'}



plt.clf()
PLT_DetectorCoverage(fillbetween=True, ax=None)
PLT_IntegrationRange(redshift=z, lower=10., upper=1E7, fillbetween=False, ax=None)

for model,version, det in zip(models, versions, dets):
	lsty = '--' if 'G' in det else '-'
	lab = '%s (v2)'%model if '2' in version else '%s'%model
	PLT_NuFnu(model, nufnutype='combined', ylims=yLims, ax=None, 
				**dict(color=colordict[model], ls=lsty, lw=1.5, alpha=0.5, label=lab)) #'%s'%model))
plt.legend(loc=8, numpoints=1, labelspacing=0.1, handletextpad=0.2)
plt.xlabel('Energy (keV)')
plt.ylabel(r'$\nu F_{\nu}$')
plt.tight_layout(pad=1)
outfilename = '/Users/KimiZ/GRBs2/python_modeling/%s_nufnu_severalmodels.pdf'%burst
plt.savefig(outfilename)
os.system('open %s'%outfilename)

#plt.show()






# 'best'            0
# 'upper right'     1
# 'upper left'      2
# 'lower left'      3
# 'lower right'     4
# 'right'           5
# 'center left'     6
# 'center right'    7
# 'lower center'    8
# 'upper center'    9
# 'center'          10

























'''






burst 		= 'bn090902462'
z 			= 1.82
model 		= 'cutoffpl+bbody+lpow' #,'cutoffpl+lpow']
versions 	= '-01-'
det 		= 'L'
detdir 		= ('GBM' if 'G' in det else 'GBMwLAT')
filename 	= '/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/%s/%s/xspec_modelreload_%s_%s_%s_.json'%(burst, detdir, model, model, version, det)
p 			= json.load(open(filename))
PLT_NuFnu(model=model, z=z, ylims=[1E0, 10**7], ax=None, **p)



model 		= 'cutoffpl+lpow'
versions 	= '-01-'
filename 	= '/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/%s/%s/xspec_modelreload_%s_%s_%s_.json'%(burst, detdir, model, model, version, det)
p 			= json.load(open(filename))
PLT_NuFnu(model=model, z=z, ylims=[1E0, 10**7], ax=None, **p)

plt.show()





burst = ['bn090902462', 'bn090902462']
z = [1.82, 1.82]
model = ['cutoffpl+bbody+lpow','cutoffpl+lpow']
versions = ['-01-','-01-']
det = 'L'
detdir = ('GBM' if 'G' in det else 'GBMwLAT')
filename = '/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/%s/%s/xspec_modelreload_%s_%s_%s_.json'%(burst, detdir, model, model, version, det)
p = json.load(open(filename))


if type(model) is list:
	print('a')


x 		= np.logspace(0.01, 8, 4000)
ylims 	= [1E0, 10**7]

modelpts = []
n=0
for mod in model.split('+'):
	modelpts.append(x * x * model_component(mod, **p)(x))
	n=n+1

#plt.clf()
plt.figure(figsize=(7,4))
plt.grid(which='major', axis='both', alpha=0.2)
plt.fill_between([.001, 10./(1.+z)], 10**(-14), 10**14, 
              color='grey', alpha=0.2)
plt.fill_between([10000./(1.+z), 10**8], 10**(-14), 
              10**14, color='grey', alpha=0.2)
for i in modelpts:
	plt.plot(x, i, color='grey', ls='--', lw=1, alpha=0.7)
plt.plot(x, sum(modelpts), color='blue', ls='-', lw=2, alpha=0.5)
plt.hlines(890, 8, 38000, color='pink', lw=6, alpha=0.3)        		# GBM
plt.hlines(600, 2E4, 3E8, color='goldenrod', lw=6, alpha=0.3)        	# LAT            
plt.xlim(1.0, 10**8)
plt.ylim(ylims[0], ylims[1])
plt.loglog()
plt.xlabel('Energy (keV)', fontsize=12)
plt.ylabel(r'$\nu F_{\nu}$', fontsize=16)
plt.show()



x = np.logspace(0.01, 8, 4000)

ylims = [1E0, 10**7]
#ylims = [None, None]

modelpts = []
n=0
for mod in model.split('+'):
	modelpts.append(x * x * model_component(mod, **p)(x))
	n=n+1

#plt.clf()
plt.figure(figsize=(7,4))
plt.grid(which='major', axis='both', alpha=0.2)
plt.fill_between([.001, 10./(1.+z)], 10**(-14), 10**14, 
              color='grey', alpha=0.2)
plt.fill_between([10000./(1.+z), 10**8], 10**(-14), 
              10**14, color='grey', alpha=0.2)
for i in modelpts:
	plt.plot(x, i, color='grey', ls='--', lw=1, alpha=0.7)
plt.plot(x, sum(modelpts), color='blue', ls='-', lw=2, alpha=0.5)
plt.hlines(890, 8, 38000, color='pink', lw=6, alpha=0.3)        		# GBM
plt.hlines(600, 2E4, 3E8, color='goldenrod', lw=6, alpha=0.3)        	# LAT            
plt.xlim(1.0, 10**8)
plt.ylim(ylims[0], ylims[1])
plt.loglog()
plt.xlabel('Energy (keV)', fontsize=12)
plt.ylabel(r'$\nu F_{\nu}$', fontsize=16)
plt.show()





def PLT_full_NuFnu(model, z, ylims, ax=None, **parameters):
	if ax is None:
		ax = plt.gca()
	if ylims is None:
		ylims 	= [None, None]

	x 		= np.logspace(0.01, 8, 4000)
	
	modelpts = []
	n=0
	for mod in model.split('+'):
		modelpts.append(x * x * model_component(mod, **parameters)(x))
		n=n+1

	PLT = ax.plot(x, sum(modelpts), color='blue', ls='-', lw=2, alpha=0.5)
	for i in modelpts:
		ax.plot(x, i, color='grey', ls='--', lw=1, alpha=0.7)

	# #plt.clf()
	# plt.figure(figsize=(7,4))
	#ax.grid(which='major', axis='both', alpha=0.2)
	ax.fill_between([.001, 10./(1.+z)], 10**(-14), 10**14, 
		color='grey', alpha=0.2)
	ax.fill_between([10000./(1.+z), 10**8], 10**(-14), 
		10**14, color='grey', alpha=0.2)
	
	ax.hlines(890, 8, 38000, color='pink', lw=6, alpha=0.3)        		# GBM
	ax.hlines(600, 2E4, 3E8, color='goldenrod', lw=6, alpha=0.3)        	# LAT            
	ax.set_xlim(1.0, 10**8)
	ax.set_ylim(ylims[0], ylims[1])
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('Energy (keV)', fontsize=12)
	ax.set_ylabel(r'$\nu F_{\nu}$', fontsize=16)
	return PLT

	#plt.show()


'''
