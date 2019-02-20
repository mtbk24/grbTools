from __future__ import division
import numpy


def percentage_error(actual, estimated):
    '''
    percentage_error(actual, estimated)
    Measures the percentage of error between the true value and an estimate of that value. Both must be measurements of the same thing and one value must be the "correct" value you are comparing a second result to. Otherwise, use percentage_difference.
    Percent Error: 
    Applied when comparing an experimental quantity, E, with a theoretical quantity, T,
    which is considered the correct value. The percent error is the absolute value of the
    difference divided by the correct value times 100.
    
        difference = (abs(estimated - actual)/(actual)) * 100.0
    
    '''
    difference = (abs(estimated - actual)/(actual)) * 100.0
    return difference



def percentage_difference(valueA, valueB):
    '''
    percentage_difference(valueA, valueB)
    Measures the difference between two values as a percentage.
    Percent Difference: Applied when comparing two experimental quantities, valueA and valueB, 
    neither of which can be considered the correct value. The percent difference is the 
    absolute value of the difference over the mean times 100.
    
        difference = (abs(valueB - valueA) / (0.5*(valueB + valueA)) ) * 100.0
    '''
    difference = (abs(valueB - valueA) / (0.5*(valueB + valueA)) ) * 100.0
    return difference


    
def residuals(ydata, ymodel):
    '''
    residuals are the: data - model
    These are unweighted residuals. The sigma residuals are weighted: (data-model)/sigma
    We don't do that here. 
    resids = [(i-j) for i,j in zip(ydata-ymodel)]

    '''
    resids = [(i-j) for i,j in zip(ydata-ymodel)]
    return resids
    


def rms(x):
    return numpy.sqrt(sum([i**2 for i in x]))/numpy.sqrt(len(x))


rms2 = lambda x: numpy.sqrt(numpy.mean(x**2))



