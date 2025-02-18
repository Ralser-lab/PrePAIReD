"""
Input: wide-format data
Output: Normalized wide-format data
"""

# Imports
import time #used for logging
import warnings
import numpy as np
import pandas as pd

from scipy.optimize import least_squares, minimize
from itertools import combinations
from sklearn.preprocessing import QuantileTransformer

# Global Variables
RANDOM_STATE = 42 # Import from global wrapper!

####################################################################
# Normalization functions
####################################################################  

## Quantile Normalization
def quantile_transform(df):
    """ 
    Apply quantile normalization to a DataFrame, transforming sample intensities
    to follow a Gaussian distribution with the same grand mean and sd
    as the original data.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        A DataFrame with rows representing samples and columns representing 
        protein features (or other measurements). The function assumes that the 
        data consists of intensity values that need to be transformed to improve 
        distribution consistency.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame with the same shape as the input, but with quantile-
        transformed intensities that follow a Gaussian distribution with the 
        same mean and standard deviation as the original data.

    Author: SB
    """

    # Apply quantile normalization on rows (samples)
    tf = QuantileTransformer(
        n_quantiles=len(df.T), output_distribution='normal', 
        random_state=RANDOM_STATE)
    df_qtransformed = pd.DataFrame(
        tf.fit_transform(df.T).T, columns=df.columns, index=df.index)

    # Adjust the scale and location to match original data's mean and sd
    original_mean = df.mean().mean()
    original_std = df.stack().std()
    df_adj = (df_qtransformed * original_std) + original_mean

    # Return the adjusted quantile-transformed data
    return df_adj

def quantile_normalize(df):
    """ 
    Applyquantile normalization (rank ordered quantile normalization) 

    Parameters:
    -----------
    df : pandas.DataFrame
        A DataFrame with rows representing samples and columns representing 
        features (or other measurements). The function assumes that the 
        data consists of intensity values that need to be transformed to improve 
        distribution consistency.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame with the same shape as the input, but with quantile-
        transformed intensities that follow a Gaussian distribution with the 
        same mean and standard deviation as the original data.

    Author: SB
    """

    # Rank the values 
    df_rank = df.rank(method="average")
    
    # Compute the average ranked value across all samples
    df_sorted_mean = df.stack().groupby(df_rank.stack().astype(int)).mean()
    
    # Replace the ranked values with the average ranked values
    df_normalized = df_rank.stack().astype(int).map(df_sorted_mean).unstack()
    
    # Return the normalized dataframe
    return df_normalized


## Median Scaling    
def median_scaling_normalization(df, axis:int=0):
    """
    Perform median scaling normalization.
    
    author: BZ
    
    Parameters
    ----------
    intensities : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples)

    Returns
    -------
    np.array
        normalization factors for all samples
    """

    # Throw a ValueError if axis is not 1 or 0
    if axis not in [0, 1]:
        raise ValueError("Axis must be either 0 (rows) or 1 (columns)")

    return np.nanmedian(np.nanmedian(df, axis=axis))/np.nanmedian(df,axis=axis)


def median_scaling_intensity(df, axis=0):
    """
    Calculate the median of median intensities, used for normalizing a test 
    set to a training set.

    Author: BZ

    Parameters
    ----------
    df : np.array or pd.DataFrame
        Matrix of intensities with rows as precursors and columns as samples.
        
    axis : int, optional (default=0)
        Axis along which to calculate the first median. Set to 0 for precursors 
        (rows) and 1 for samples (columns).

    Returns
    -------
    float
        Median of the median intensities.

    """
    
    if axis not in [0, 1]:
        raise ValueError("Axis must be either 0 (rows) or 1 (columns)")

    return np.nanmedian(np.nanmedian(df, axis=axis))

def median_normalization(df):
    """
    median normalization 
    substracting median of every column, resulting in column medians of 0.
    
    author: BZ
    
    Parameters
    ----------
    df : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples)

    Returns
    -------
    intensities : np.array
        matrix of normalized intensities with column medians of 0 
        (1st dim peptides, 2nd dim samples)

    """
    return df.sub(df.median(axis=0),axis=1)


def mean_normalization(df):
    """
    mean normalization 
    (substracting mean of every column, resulting in column means of 0).
    
    author: BZ
    
    Parameters
    ----------
    df : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples)

    Returns
    -------
    intensities : np.array
        matrix of normalized intensities with column means of 0 
        (1st dim peptides, 2nd dim samples)

    """
    return df.sub(df.mean(axis=0),axis=1)



####################################################################
# functions for delayed normalization (Cox et al., DOI: 10.1074/mcp.M113.031591)
# require more testing, not tested sufficiently for fractionated data, 
# poorly compared to other normalization strategies
####################################################################

def logfcsq(x, pept, df, a, b):
    """
    compute squared log fold changes as inner part of delayed normalization
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    pept : int
        index in intensities for given peptide.
    df : np.array
        matrix of intensities 
        (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    a : int
        sample a.
    b : int
        sample b.
    """
    warnings.warn(
        "This is preliminary code and hasn't been sufficiently tested.", 
        UserWarning)
    return abs(
            np.log2(
                np.nansum(
                    x[a]*df[pept, a, :])/np.nansum(x[b]*df[pept, b, :])))**2


 
def Hp(x, pept, df):
    """
    compute H_peptide(x) as middle part of delayed normalization
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    pept : int
        index in intensities for given peptide.
    df : np.array
        matrix of intensities 
        (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    """
    warnings.warn(
        "This is preliminary code and hasn't been sufficiently tested.", 
        UserWarning)
    
    return np.nansum(
            list(
                map(
                    lambda xx: 0 if df[
                        pept, xx[0], :].sum() == 0 or df[pept, xx[1], :].sum() == 0 else logfcsq(
                        x, pept, df, xx[0], xx[1]), 
                        combinations(range(len(x)), 2))))


def H(x, df, pept):
    """
    ### ADDED PEPT AS INPUT, CD
    compute H(x) for all peptides
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    df : np.array
        matrix of intensities 
        (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    """
    warnings.warn(
        "This is preliminary code and hasn't been sufficiently tested.", 
        UserWarning)
    return [Hp(
            x, pept, df) for peptide in range(
                df.shape[0])]


def Hsum(x, df, pept):
    """
    ### ADDED PEPT AS INPUT, CD

    calculate sum of H functions for every peptide.
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    df : np.array
        matrix of intensities 
        (1st dim peptides, 2nd dim samples, 3rd dim fractions).


    """
    return np.sum(
        [Hp(x, pept, df) for peptide in range(df.shape[0])])


# TODO
# logging inside this function, because more information to be logged than for 
# other normalization strategies. 
# Could be moved outside function, but would change output.
# def delayed_normalization(
#         df, method="trf", bounds=(0,np.inf), logfilename=None):
#     """
#     calculate delayed normalization coefficients for given matrix of intensities 
#     as introduced by Cox et al., DOI: 10.1074/mcp.M113.031591 
    
#     author: BZ
    
#     Parameters
#     ----------
#     df : np.array
#         matrix of intensities 
#         (1st dim peptides, 2nd dim samples, 3rd dim fractions).
#     method : str, optional
#         Algorithm for least squares minimization 
#         (see scipy.optimize.least_squares), 
#         one of {‘trf’, ‘dogbox’, ‘lm’}. The default is "trf".
#     bounds : tuple, optional
#         Upper and lower bound for minimization algorithm. 
#         The default is (0,np.inf).

#     Returns
#     -------
#     res.x : np.array
#         normalization factors for all samples and fractions (if present).

#     """
#     warnings.warn(
#         "This is preliminary code and hasn't been sufficiently tested.", 
#         UserWarning)

#     #global intensities
#     x0 = np.ones(df.shape[1:]).flatten()
#     if logfilename is not None:
#         logfile=open(logfilename,"a")
#         logfile.write(
#             time.strftime(
#                 "%Y-%m-%d%t%H:%M:%S")+"\tBegin delayed normalization...\n")
#         logfile.close()
#     if method in ["trf","dogbox","lm"]:
#         res = least_squares(
#             fun=H, x0=x0, method=method, bounds=bounds, xtol=0.01, args=[df])  
#         # Trust-Region-Reflective Algorithm
#     else:
#         #bounds=[(0.1,1)]*len(x0)
#         bounds=[bounds]*len(x0)
#         res=minimize(
#             fun=Hsum, x0=x0, method=method, bounds=bounds, 
#             options={"disp":True, "maxiter":100*df.shape[0]}, args=df)
#     res = least_squares(
#         fun=H, x0=x0, method=method, bounds=(0, np.inf), xtol=0.01, args=[df])  
#     # Trust-Region-Reflective Algorithm
#     res = least_squares(
#         fun=H,x0=x0,method="lm",max_nfev=2,ftol=10000,args=[df]) 
#     #Levenberg-Marquardt
#     if logfilename is not None:
#         logfile=open(logfilename,"a")
#         logfile.write(
#             time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tEnd normalization\n")
#         logfile.write(
#             "\t\t\t\tReason for minimization termination: "+res.message+"\n")
#         try:
#             logfile.write("\t\t\t\t"+res.nfev + "\n")
#         except:
#             pass
#         logfile.close()
#     return(res.x)


####################################################################

####################################################################
#logging EXAMPLE in main
####################################################################

# def main(intensities,logfilename=None):
#     if logfilename is not None:
#         logfile=open(logfilename,"a")
#         logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tNormalizing using median scaling normalization.\n")
#         logfile.close()
#     median_scaling_intensity(intensities,axis=0)
#     if logfilename is not None:
#         logfile=open(logfilename,"a")
#         logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tFinished normalization.\n")
#         logfile.close()
    
def main():
    pass