"""
Input: wide-format data
Output: Normalized wide-format data
"""
import numpy as np
import pandas as pd

def example():
    pass

def logtransform(df, base = 2):
    # check for missing vals
    if df.isna().any().any() > 0:
        # abort
        print('missing values')
        return df
    
    # check for 0 or negative values
    if (df < = 0).any().any() > 0:
        # abort
        print('negative and/or zero values')
        return df
    
    # return log transformed df
    return np.log(df)/np.log(base)

def quantile_norm():
    pass 

def z_norm(df): 
    #calculate mean
    dfmean = df.mean(skipna = True)
    #calculate variance
    dfvar = df.var(skipna = True)
    #standardize df
    df_z = (df - dfmean)/np.sqrt(dfvar)
    #return z-normalized data
    return df_z 

def main():
    return example()
