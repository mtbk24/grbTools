from __future__ import division
import numpy as np
from scipy.optimize import minimize



def LogLikelihood(parameters, x, y, stat='gaussian'):
    '''
    LogLikelihood(parameters, x, y, stat='gaussian')
    
    uses: from scipy.optimize import minimize
    
    This is currently only a Gaussian Log Likelihood Function.
    
    LogLikelihood(parameters, x, y)
    
    parameters: m,b,sigma and returned in that order.
    
    result = minimize(LogLikelihood, np.array([1,1,1]), method='L-BFGS-B', args=(x, y)); result
    
    result = minimize(LogLikelihood,
                  np.array([1,1,0.2]), 
                  method    = 'L-BFGS-B',
                  args      = (x, y))
                  
    '''
    m, b, sigma       = parameters
    ymodel            = m * x + b
    n                 = float(len(ymodel))
    error             = y - ymodel  # y data - linear model
    L   = ((n/2.) * np.log(2.*np.pi*sigma**2) + 1./(2.*sigma**2)* np.dot(error.T,error))
    return L


def Prior(parameters):
    '''
    Prior(parameters)
    '''
    m, b, f = parameters
    if b < 0:
        if 0.0 < m < 1.0 and -40 < b < -10 and -1.0 < f < 1.0:
            return 0.0
    else:
        if -1.0 < m < 1.0 and 0 < b < 10.0 and -1.0 < f < 1.0: #0.0 < f < 1.0:
            return 0.0
    return -np.inf


def Prob(parameters, x, y):
    '''
    Prob(parameters, x, y)
    '''
    lp = Prior(parameters)
    if not np.isfinite(lp):
        return -np.inf
    return lp - LogLikelihood(parameters, x, y)





