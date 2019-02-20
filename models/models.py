from __future__ import division
import numpy as np
from numpy import exp, power, float64, array, inf, logical_and, log, log10, exp
from math import atanh
from scipy.integrate import quad, quadrature


def BAND_RMFIT( energy, alpha, beta, epeak, norm):
    '''
    BAND_RMFIT( energy, alpha, beta, epeak, norm)

    RMFIT Band Function.  Takes Epeak values only. Can not use with XSPEC results
    unless the 'tem' was converted to an 'epeak' value. 

    This is the same function as XSPEC's, but with the 'epeak' replacing 'tem'.

    tem = epeak/(2.+alpha) 
    
    '''
    N = norm
    Ep = epeak
    alpha = alpha
    beta = beta
    eng = energy

    cond1 = eng < (alpha-beta)*Ep/(2.+alpha)
    cond2 = eng >= (alpha-beta)*Ep/(2.+alpha)
    
    band = np.piecewise(eng, [cond1, cond2],\
            [lambda x: N * ( power(x/100., alpha) * exp(-x*(2.+alpha)/Ep) ), \
             lambda x: N * ( power( (alpha-beta)*Ep/(100.*(2.+alpha)), alpha-beta) * exp(beta-alpha)*power(x/100., beta))])

    return band


def BAND(energy, alpha, beta, enterm, norm, entype='E0'):
    '''
    XSPEC Band Function.
    
    entype defines what you are passing to the 'enterm'.
    All enterm's will be converted into E0, which is represented by the term t.
    This is 'tem' in the XSPEC funciton, aka, the high-energy cutoff.
    
    The relations with ebreak and epeak:
    ebreak = ((alpha - beta) * epeak) / (alpha + 2.)
    
    and epeak and E0:
    E0 = epeak / (alpha + 2.)
    
    
    if entype = 'epeak', you are passing epeak and is converted to E0, which 
        is the t parameter in the function.
    if entype = 'ebreak', you are passing ebreak and it is converted into E0.
    if entype = 'E0', nothing happens because function is set to take E0.
    
    '''
    eng     = energy
    N       = norm
    a       = alpha
    b       = beta
    
    if entype == 'epeak':
        t   = enterm / (2.+a)  # E0 = epeak / (2.+alpha)
    elif entype == 'E0':
        t   = enterm           # E0 = E0
    elif entype == 'ebreak':
        t   = enterm/(a-b)  # E0 = ebreak / (alpha - beta)
    else:
        raise Exception, "entype must be 'epeak', 'E0', or 'ebreak' "

    cond1 = eng < (a-b) * t
    cond2 = eng >= (a-b) * t
    
    band = np.piecewise(eng, [cond1, cond2],\
            [lambda x: N * (power(x/100., a) * exp(-x/t)), \
             lambda x: N * (((((a-b)*t)/100.0)**(a-b)) * (exp(b-a))*((x/100.0)**b))])

             #lambda x: N * (power((a-b)*Ep/(100.*(2.+alpha)), a-b) * exp(b-a) * power(eng/100.,b))])

    return band
#
#    if eng < ( (a-b) * t ):
#        return  N * (((eng/100.0)**a) * (exp(-eng/t)))
#    else:
#        return  N * (((((a-b)*t)/100.0)**(a-b)) * (exp(b-a))*((eng/100.0)**b))


def BBODY_RMFIT(energy, kT, norm):
    '''
    BBODY_RMFIT(energy, kT, norm)

    GBM (RMFIT?) Blackbody

    eqn = N * power(eng,2) * power(exp(eng/temp)-1, -1)

    '''
    eng     = energy
    N       = norm
    temp    = kT
    return N * power(eng,2) * power(exp(eng/temp)-1, -1)


def BBODY(energy, kT, norm):
    '''
    XSPEC Blackbody
    '''
    eng     = energy
    temp    = kT
    N       = norm
    return N * (((eng**2)*(8.0525)) / ((temp**4) * (exp(eng/temp)-1)))
    #if eng <= (709.666 * kT): # to avoid exp overflow error
    #    return N * (((eng**2)*(8.0525)) / ((kT**4) * (exp(eng/kT)-1)))
    #else:
    #    return 0


def LPOW_RMFIT(energy, index, norm):
    '''
    LPOW_RMFIT(energy, index, norm)

    RMFIT power-law.

    Should be same as XSPEC's PL model.
    
    eqn = A*(x/Epiv)**index
    '''
    N      = norm
    eng    = energy
    return N * ((eng/100.0)**index)


def LPOW(energy, index, norm):
    '''
    XSPEC power-law.  
    Should be same as the RMFIT PowerLaw.
    '''
    N      = norm
    eng    = energy
    return N * ((eng/100.0)**index)


def COPL_RMFIT(energy, index, epeak, norm):
    '''
    COPL_RMFIT(energy, index, epeak, norm)

    RMFIT's copl model.
    '''
    eng     = energy
    N       = norm
    Ep      = epeak
    return N * exp(-eng * (2.0+index)/Ep ) * power(eng/100.0, index)


def COPL(energy, index, enterm, norm, entype='E0'):
    '''
    XSPEC Cutoff Power-Law.

    Takes only 'E0' or 'epeak' as entype. 
    E0 is also known as:   high-energy cutoff, e-folding energy, tem. 

    
    POSSIBLE PROBLEMS WITH THIS FUNCTION AND XSPEC PARAMETER VALUES:

        WE ADJUSTED THE PARAMETERS FOR 'cutoffpl' TO BE USED WITH THE
        FUNCTION NORMALIZED TO 100 KEV (THIS ONE).  SINCE THE XSPEC FUNCTION
        IS NOT NORMALIZED TO 100, BUT INSTEAD TO 1.,  WE MUST RE-ADJUST THE 
        PARAMETERS SO THEY CAN BE USED WITH THIS FUNCTION.
        ONCE YOU
        FLIP SIGN ON PHOTON INDEX AND MULTIPLY THE NORMALIZATION BY 100^(-PHOTON INDEX)
        
        if 'cutoffpl' in model:
            print("\n \n !*!*!*! FIXING CUTOFFPL PARAMS !*!*!*! \n \n ")
            for i,name in enumerate(parNames):
                if 'PhoIndex__1' in name:
                    globals()[name]     = -pyxResults[name][0]
                elif 'norm__3' in name:
                    globals()[name]     = (pyxResults[name][0])*(100.**(-pyxResults['PhoIndex__1'][0]))
                else:
                    globals()[name]     = pyxResults[name][0]
        else:
            for i,name in enumerate(parNames):
                globals()[name]     = pyxResults[name][0]
                
    
    '''
    
    if entype == 'epeak':
        t = enterm / (2.+index)  # E0 = epeak / (2.+alpha)
    elif entype == 'E0':
        t = enterm
    else:
        raise Exception, "entype must be either 'E0' or 'epeak' "
    N = norm
    eng = energy
    return N * ((eng/100.0)**index) * (exp(-eng/t))


def SBPL(energy, alpha, beta, enterm, norm, entype='ebreak'):
    '''
    entype defines what you are passing to the 'enterm'.
    
    if entype = 'epeak', you are passing epeak and ebreak must be found to 
        run this function.  
    if entype = 'ebreak', then conversion is skipped because this function is 
        setup to take ebreak.

    
    Conversion between epeak and ebreak:
    # epeak = (ebreak) * (10**(0.3 * atanh((a+b+4)/(a-b))))
    # ebreak = (epeak) / (10**(0.3 * atanh((a+b+4)/(a-b))))
    # where a and b are alpha and beta
    '''

    eng  = energy
    a    = alpha
    b    = beta
    d    = 0.3
    N    = norm
    
    if entype == 'epeak':
        k = (enterm) / (10**(d * atanh((a+b+4.)/(a-b))))
        # epeak = (ebreak) * (10**(0.3 * atanh((a+b+4)/(a-b))))
        # ebreak = (epeak) / (10**(0.3 * atanh((a+b+4)/(a-b))))
    elif entype == 'ebreak':
        k = enterm  # THIS FUNCTION USES EBREAK
    else:
        raise Exception, "entype must be either 'epeak' or 'ebreak' "
    p1   = (beta - alpha)/2.
    p2   = (alpha + beta)/2.
    p3   = ( log10(100.0/k)/d )
    p4   = ( log10(eng/k)/d )
    return N * ((eng/100.0)**p2) * (10**((p1 * d * log((exp(p4) + exp(-p4))/2.)) - (p1 * d * log((exp(p3) + exp(-p3))/2.))))


def SBPL_RMFIT(energy, alpha, beta, ebreak, norm):
    '''
    SBPL_RMFIT(energy, alpha, beta, ebreak, norm)

    RMFIT's SBPL model. Should be same as XSPEC's, but XSPEC's 
    has more options for the characteristic energy term.


    Conversion between epeak and ebreak:
    epeak = (ebreak) * (10**(0.3 * atanh((a+b+4)/(a-b))))
    ebreak = (epeak) / (10**(0.3 * atanh((a+b+4)/(a-b))))
           where a and b are alpha and beta

    '''

    eng     = energy
    alpha   = alpha
    beta    = beta
    ebkscl  = 0.3
    N       = norm
    Ebk     = ebreak
    
    p1   = (beta - alpha)/2.
    p2   = (alpha + beta)/2.
    p3   = ( log10(100.0/Ebk)/ebkscl )
    p4   = ( log10(eng/Ebk)/ebkscl )
    return N * ((eng/100.0)**p2) * (10**((p1 * ebkscl * log((exp(p4) + exp(-p4))/2.)) - (p1 * ebkscl * log((exp(p3) + exp(-p3))/2.))))



