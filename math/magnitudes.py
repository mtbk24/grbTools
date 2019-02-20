from __future__ import division

def magnitude(x):
    '''
    magnitude(x)
    
    x - float of a number.
    
    This function returns the order of magnitude of x.
    This is a good tool for making Latex Tables.
    
    x = 1.92795239E-9
    magnitude(x) returns 1e-09
    now,   x/magnitude(x)   will be 1.9279523899999997
    
    x = 34.295E+5
    magnitude(x) returns  1000000
    now,   x/magnitude(x)   will be 3.4295
    
    '''
    import numpy
    if x > 1.0:
        return 1*10**( int(numpy.log10(x)) )
    else:
        return 1*10**( int(numpy.log10(x)) - 1)

