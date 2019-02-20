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
import numpy as np
from scipy import integrate
from math import pi
from numpy import exp, power, float64, array, inf, logical_and, log, log10, exp
from math import atanh
from scipy.integrate import quad, quadrature
import json




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
			
	
	'''
	a 			= parsdict['PhoIndex'][0]
	t 			= parsdict['HighECut'][0]
	N 			= parsdict['norm'][0]
	eng 		= energy
	return N * ((eng/100.0)**a) * (exp(-eng/t))














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

   









# burst = 'bn'

# z = 2.83

# model = 'sbpl+bbody'

# version = '-01-'

# det = 'L'


# bursts = ['bn080916009', 'bn090323002', 'bn090328401', 'bn090510016', 
# 	'bn090902462', 'bn090926181', 'bn091003191', 'bn091208410', 
# 	'bn100728095', 'bn110731465', 'bn130518580', 'bn131108862', 'bn131231198']

# zs = [4.35, 3.57, 0.736, 0.903, 1.822, 2.1062, 0.8969, 1.0633, 1.567, 2.83, 2.49, 2.4, 0.642]


# bursts = ['bn131108862',]

# zs = [2.4,]


# models = ['grbm+lpow',]*len(bursts)

# versions = ['-01-',]*len(bursts)

# det = 'L'



bursts = ['bn090902462','bn090902462']

zs = [1.82, 1.82]


models = ['cutoffpl+bbody+lpow','cutoffpl+lpow',]

versions = ['-01-',]*len(bursts)

det = 'L'



# models = ['sbpl+lpow','grbm+lpow','sbpl+bbody+lpow','grbm+bbody+lpow','cutoffpl+lpow','cutoffpl+bbody+lpow']

# bursts = ['bn090510016',]*len(models)

# zs = [0.903,]*len(models)

# versions = ['-01-',]*len(models)

# det = 'L'


for burst,z,model,version in zip(bursts,zs,models,versions):
	detdir = ('GBM' if 'G' in det else 'GBMwLAT')
	filename = '/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/%s/%s/xspec_modelreload_%s_%s_%s_.json'%(burst, detdir, model, model, version, det)
	p = json.load(open(filename))
	fluxes = []
	for mod in model.split('+'):
		# REST-FRAME
		flux = Calc_Flux_FromFile(emin=10, emax=10000, redshift=z, model=mod, **p)
		#
		# OBSERVER FRAME
		#flux = Calc_Flux_FromFile(emin=.1, emax=100, redshift=None, model=mod, **p)
		fluxes.append(flux)
	fullflux = sum(fluxes)
	print('\n')
	print(p)
	print('%s  %f'%(burst, z))
	print('%s  %s'%(model, version))
	print('')
	print(fluxes)
	print(fullflux)
	print('')
	for i,j in enumerate(model.split('+')):
		print('%8.2f  %s'%( ((fluxes[i]/fullflux)*100), j))
	print('\n')
	print('--'*30)








for burst,z,model,version in zip(bursts,zs,models,versions):
	detdir = ('GBM' if 'G' in det else 'GBMwLAT')
	filename = '/Users/KimiZ/GRBs2/analysis/LAT/%s/PYXSPEC/%s/%s/xspec_modelreload_%s_%s_%s_.json'%(burst, detdir, model, model, version, det)
	p = json.load(open(filename))
	fluxes = []
	for mod in model.split('+'):
		# REST-FRAME
		flux = Calc_Flux_FromFile(emin=10, emax=10000, redshift=z, model=mod, **p)
		#
		# OBSERVER FRAME
		#flux = Calc_Flux_FromFile(emin=.1, emax=100, redshift=None, model=mod, **p)
		fluxes.append(flux)
	fullflux = sum(fluxes)
	print('\n')
	print(p)
	print('%s  %f'%(burst, z))
	print('%s  %s'%(model, version))
	print('')
	print(fluxes)
	print(fullflux)
	print('')
	for i,j in enumerate(model.split('+')):
		print('%8.2f  %s'%( ((fluxes[i]/fullflux)*100), j))
	print('\n')
	print('--'*30)

