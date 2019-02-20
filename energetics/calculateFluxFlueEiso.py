from __future__ import division
from scipy import integrate
from math import pi
from Zoldak.Cosmology.luminositydistance import LumDist
import spectralTools
from spectralTools.models import BAND, SBPL, COPL, BBODY, LPOW, BAND_RMFIT, SBPL_RMFIT, COPL_RMFIT, BBODY_RMFIT, LPOW_RMFIT


def Calc_Flux(model, Pars, emin, emax, redshift=None):
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
    listofmodels = ['grbm', 'sbpl', 'cutoffpl', 'bbody', 'lpow']

    keVtoerg      = 1.60217657E-9
    if redshift is not None:
        emin      = emin/(1.+redshift)
        emax      = emax/(1.+redshift)
    else:
        print('\n *** WARNING: You are using observer-frame fluxes. Eiso requires rest-frame fluxes. Please use redshift. *** \n')
        emin     = emin
        emax     = emax
    
    P = Pars
    
    if model == 'grbm':
        function = lambda energy: energy * BAND(energy,  
                                                P['grbm']['alpha'], 
                                                P['grbm']['beta'], 
                                                P['grbm']['enterm'], 
                                                P['grbm']['norm'], 
                                                entype=P['grbm']['entype']
                                                )
    
    elif model == 'sbpl':
        function = lambda energy: energy * SBPL(energy,  
                                    P['sbpl']['alpha'], 
                                    P['sbpl']['beta'], 
                                    P['sbpl']['enterm'], 
                                    P['sbpl']['norm'], 
                                    entype=P['sbpl']['entype'])
    
    elif model == 'cutoffpl':
        function = lambda energy: energy * COPL(energy,  
                                    P['cutoffpl']['alpha'], 
                                    P['cutoffpl']['enterm'], 
                                    P['cutoffpl']['norm'], 
                                    entype=P['cutoffpl']['entype'])
                                                
                                                
    
    elif model == 'lpow':
        function = lambda energy: energy * LPOW(energy,  
                                    P['lpow']['alpha'], 
                                    P['lpow']['norm'])
                                                
    
    
    elif model == 'bbody':
        function = lambda energy: energy * BBODY(energy,  
                        P['bbody']['kT'], 
                        P['bbody']['norm'])
    else:
        raise Exception, "Don't recognize models. Options are: %r"%listofmodels

    nrgFlux = integrate.quad(function, emin, emax, limit=100)[0] * keVtoerg
    return nrgFlux


def Calc_Fluence(model, Pars, emin, emax, duration, redshift=None):
    if redshift:
        flux        = Calc_Flux(model, Pars, emin, emax, redshift)
    else:
        flux        = Calc_Flux(model, Pars, emin, emax)
    fluence     = flux * duration
    return fluence



def Calc_Eiso(model, Pars, emin, emax, duration, redshift):
    flux        = Calc_Flux(model, Pars, emin, emax, redshift)
    fluence     = flux * duration
    dL          = LumDist(redshift)
    eiso        = ((4.0 * pi * (dL**2))/(1.+redshift)) * fluence
    return eiso






# def Calc_Flux(model, Pars, emin=10., emax=10000., redshift=None):
#     '''
#     Calc_Flux(model, Pars, emin=10., emax=10000., redshift=None)
    
#     model:  str, model name.
#             ex:  'grbm', 'grbm+bbody', 'grbm+bbody+lpow'
#             all individual model options:
#                 grbm, sbpl, cutoffpl, lpow, and bbody
    
#     Pars:   dict, model parameters.  SEE BELOW.
    
#     emin:   float, lower energy integration limit. Do not k-correct yourself.
    
#     emax:   float, upper energy integration limit. Do not k-correct yourself.
    
#     redshift: float, redshift of GRB for k-correcting. If not k-correcting, do not enter anything.
    
#     kcorrect:   [True|False]  
#                 True will use k-correction: emin/(1.+redshift) to emax/(1.+redshift)
    
    
    
#     # REGARDLESS OF WHETHER THERE IS ONE MODEL OR MORE THAN ONE, PARAMETER DICTIONARIES SHOULD
#     # ALWAYS BE SET UP THIS WAY.  THIS IS TO ENSURE PARAMETERS BETWEEN TWO DIFFERENT MODELS DO NOT GET MIXED UP.
    
#     pars = {'grbm':     {'alpha':   -1.0,
#                          'beta':    -2.5,
#                          'enterm':  300.0,
#                          'norm':    0.01,
#                          'entype':  'epeak'},
#             'sbpl':     {'alpha':   -1.0,
#                          'beta':    -2.5,
#                          'enterm':  300.0,
#                          'norm':    0.01,
#                          'entype':  'epeak'},
#             'cutoffpl': {'alpha':   -1.0,
#                          'enterm':  300.0,
#                          'norm':    0.01,
#                          'entype':  'epeak'},
#             'lpow':     {'index':   -1.0,
#                          'norm':    0.001},
#             'bbody':    {'kT':      10.0,
#                          'norm':    0.001} } 
                     
    
#     '''
#     keVtoerg    = 1.60217657E-9
#     if redshift is not None:
#         emin_K      = emin/(1.+redshift)
#         emax_K      = emax/(1.+redshift)
#     else:
#         print('\n *** WARNING:  NOT USING K-CORRECTED ENERGIES. *** \n')
#         emin_K      = emin
#         emax_K      = emax
    
#     models = model.split('+')
#     P = Pars
    
    
#     nrgFlux = []
#     for mod in models:
#         if 'grbm' in mod:
#             function = lambda energy: energy * BAND(energy,  
#                                                     P['grbm']['alpha'], 
#                                                     P['grbm']['beta'], 
#                                                     P['grbm']['enterm'], 
#                                                     P['grbm']['norm'], 
#                                                     entype=P['grbm']['entype'])
            
#             #integrate.quad(function, emin_K, emax_K)[0] * keVtoerg

        
#         if 'sbpl' in mod:
#             function = lambda energy: energy * SBPL(energy,  
#                                         P['sbpl']['alpha'], 
#                                         P['sbpl']['beta'], 
#                                         P['sbpl']['enterm'], 
#                                         P['sbpl']['norm'], 
#                                         entype=P['sbpl']['entype'])
        
#         if 'cutoffpl' in mod:
#             function = lambda energy: energy * COPL(energy,  
#                                         P['cutoffpl']['alpha'], 
#                                         P['cutoffpl']['enterm'], 
#                                         P['cutoffpl']['norm'], 
#                                         entype=P['cutoffpl']['entype'])
                                                    
                                                    
        
#         if 'lpow' in mod:
#             function = lambda energy: energy * LPOW(energy,  
#                                         P['lpow']['alpha'], 
#                                         P['lpow']['norm'])
                                                    
        
        
#         if 'bbody' in mod:
#             function = lambda energy: energy * BBODY(energy,  
#                             P['bbody']['kT'], 
#                             P['bbody']['norm'])

#         Flux = integrate.quad(function, emin_K, emax_K, limit=100)[0] * keVtoerg
#         nrgFlux.append(Flux)
#     nrgFlux      = sum(nrgFlux)
#     return nrgFlux

    
