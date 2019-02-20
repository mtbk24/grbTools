from __future__ import division
import numpy
#from scipy import polyfit, polyval

def line(x, m, b):
    '''
    line(x, m, b)
    
    y = mx + b
    m = slope
    b = y-intercept
    
    '''
    return m*x + b


def polynomial(x, y, polyorder):
    '''
    polynomial(x, y, polyorder)
    polyorder - polynomial order 0-5.
    
    result = polynomial(x, y, 1)
    
    x - array of x values to fit polynomial to.
    y - array of y values to fit polynomial to.


    fits a polynomial of degree = polyorder to the x and y data to find a best fit.
    returns the funciton of best fit to be used to generate new x and y data.


    New x and y data based on the polynomial returned as best fit:
    func = polynomial(x, y, polyorder = 1)
    x_new       = numpy.linspace(min(x), max(x), 100)
    y_new       = func(x_new)
    plt.plot(x_new, y_new, color = 'grey', ls='--', lw=3, alpha = 0.7)

    or linear fit to data in log-log space:
    xx          = numpy.linspace(-2, 2, 100)
    yy          = func(xx)
    plt.plot(10**xx, 10**yy, color = 'grey', ls='--', lw=3, alpha = 0.7)
      
    contents of function:
    results     = numpy.polyfit(x, y, polyorder)
    Func        = numpy.poly1d(results)
    return Func
        
    '''
    results     = numpy.polyfit(x, y, polyorder)
    Func        = numpy.poly1d(results)
    return Func


def exponential_function(x, a, b, c):
    '''
    exponential_function(x, a, b, c)
    
    
    x       = numpy.asarray(x)
    y       = numpy.asarray(y)
    
    # coefficients and covariance
    popt, pcov = scipy.optimize.curve_fit(exponential_function, x, y)
    
    # exponential fit to variables in log-linear space
    xdata   = numpy.linspace(-2, 2, 10000)
    ydata   = exponential_function(xdata, *popt)
    plt.plot(10**xdata, ydata, color='hotpink', lw=2.5, alpha=0.7)
    '''
    return a * numpy.exp(-b * x) + c  # GOOD



def rms(x):
    return numpy.sqrt(sum([i**2 for i in x]))/numpy.sqrt(len(x))


rms2 = lambda x: numpy.sqrt(numpy.mean(x**2))

