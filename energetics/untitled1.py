from __future__ import division
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
from scipy import integrate
from math import pi
from Zoldak.Cosmology.luminositydistance import LumDist
import spectralTools
from spectralTools.xspecmodels import BAND, SBPL, COPL, BBODY, LPOW


#def Calc_Flux(model, Pars, emin, emax, redshift=None):


def Calc_Flux(emin, emax, redshift, model, **Pars):
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
        #emin     = emin
        #emax     = emax
    if model == 'grbm':
        function = lambda energy: energy * BAND(energy, **Pars)
    
    elif model == 'sbpl':
        function = lambda energy: energy * SBPL(energy, **Pars)
    
    elif model == 'cutoffpl':
        function = lambda energy: energy * COPL(energy, **Pars)
                                                
    elif model == 'lpow':
        function = lambda energy: energy * LPOW(energy, **Pars)
    
    elif model == 'bbody':
        function = lambda energy: energy * BBODY(energy, **Pars)
    else:
        raise Exception, "Don't recognize models. Options are: %r"%listofmodels

    nrgFlux = integrate.quad(function, emin, emax, limit=100)[0] * keVtoerg
    return nrgFlux
















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
        function = lambda energy: energy * BAND(energy, **Pars)
    
    elif model == 'sbpl':
        function = lambda energy: energy * SBPL(energy, **Pars)
    
    elif model == 'cutoffpl':
        function = lambda energy: energy * COPL(energy, **Pars)
                                                
    elif model == 'lpow':
        function = lambda energy: energy * LPOW(energy, **Pars)
    
    elif model == 'bbody':
        function = lambda energy: energy * BBODY(energy, **Pars)
    else:
        raise Exception, "Don't recognize models. Options are: %r"%listofmodels

    nrgFlux = integrate.quad(function, emin, emax, limit=100)[0] * keVtoerg
    return nrgFlux
	# function = eval('lambda energy: energy * %s(energy, **Pars)'%modelfunc)
	# nrgFlux = integrate.quad(function, emin, emax, limit=100)[0] * keVtoerg
	# return nrgFlux
    



# p = {'grbm':{'alpha':-1.029585538992228, 
# 	 		'beta':-2.198223234704878, 
# 	 		'enterm': 526.2095854654398, 
# 	 		'norm': 0.01752245438164585, 
# 	 		'entype': 'E0'},
# 	 'bbody':{},
# 	 'lpow':{}}


p = {
    "grbm": {
        "alpha": -1.1316211655148538,  
        "beta": -2.2785592988736734, 
        "tem": 1108.9523058886602, 
        "norm": 0.010934933257653946}, 
    "bbody": {
        "kT": 44.539900277976834, 
        "norm": 1.7123658873772487}, 
    "lpow": {
        "alpha": -1.886795470530438, 
        "norm": 0.0009080978250198132}}

Calc_Flux(emin=10, emax=10000, redshift=4.35, model='grbm', **p)



















def Calc_Fluence(model, Pars, emin, emax, duration, redshift=None):
    flux        = Calc_Flux(model, Pars, emin, emax, redshift)
    fluence     = flux * duration
    return fluence



def Calc_Eiso(model, Pars, emin, emax, duration, redshift=None):
    #fluence = Calc_Fluence(model, Pars, emin, emax, duration, redshift)
    flux        = Calc_Flux(model, Pars, emin, emax, redshift)
    fluence     = flux * duration
    dL          = LumDist(redshift)
    eiso        = ((4 * pi * (dL**2))/(1.+redshift)) * fluence
    return eiso




'/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/GBMwLAT/%s/Plots_%s_%s_L_.ps'%(burst, model, model, version)


for burst, model, version in zip(k.trigger.tolist(), k.function.tolist(), k.version.tolist()):
	print('/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/GBMwLAT/%s/Plots_%s_%s_L_.ps'%(burst, model, model, version)
)
