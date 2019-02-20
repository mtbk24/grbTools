
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, kendalltau, spearmanr, ks_2samp


def correlation_matrix(dataframe, method):
    '''
    correlation_matrix(dataframe, method):
    -- Uses Pandas Dataframes.
    -- must remove all columns with strings and NaN.

    dataframe: pandas data frame with floats and ints (numbers) only.

    method:   'pearsonr', 'spearmanr', 'kendalltau'

    -- uses scipy's pearsonr, spearmanr, and kendalltau
    -- uses dataframe column names and axes.
    -- returns 2 dataframes, a coefficient and a p-value.
    
    coeff, pvals = correlation_matrix(dataframe=,
                               method='spearmanr')
    
    '''
    df = dataframe.copy() 
    coeffmat    = np.zeros([df.shape[1], df.shape[1]]) 
    pvalmat     = np.zeros([df.shape[1], df.shape[1]]) 
    colnames    = df.columns
    for i,col1 in enumerate(colnames):
        for j,col2 in enumerate(colnames):
            # eval performs: pearsonr(a,b) or other stat.
            # method must be pearsonr, not pearson b/c scipy's 
            # func is pearsonr.
            corrtest = eval('%s(df.eval(col1), \
                            df.eval(col2))'%method)
            coeffmat[i,j]   = corrtest[0]
            pvalmat[i,j]    = corrtest[1]
    coeffmat    = pd.DataFrame(coeffmat, 
                               columns=colnames, index=colnames)        
    pvalmat     = pd.DataFrame(pvalmat, 
                               columns=colnames, index=colnames) 
    return coeffmat, pvalmat



