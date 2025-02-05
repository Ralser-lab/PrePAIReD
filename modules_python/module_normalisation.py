"""
Input: wide-format data
Output: Normalized wide-format data
"""
import numpy as np
import pandas as pd

def example():
    pass

def logtransform():
    pass

def quantile_norm():
    pass 

def z_norm(df): 
    #calculate mean
    mean = df.mean(skipna = True)
    #calculate variance
    var = df.var(skipna = True)
    #standardize df
    df_z = (df - mean)/np.sqrt(var)
    #return z-normalized data
    return df_z 

def main():
    return example()
