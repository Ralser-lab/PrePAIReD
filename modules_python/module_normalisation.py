"""
Input: wide-format data
Output: Normalized wide-format data
"""

# Imports
import warnings
import numpy as np
import pandas as pd
import sys
import helper
import json
import os
from scipy.optimize import least_squares, minimize

from itertools import combinations
from sklearn.preprocessing import QuantileTransformer

# Global Variables
RANDOM_STATE = 42 # Import from global wrapper!

####################################################################
# Normalization functions
####################################################################  

## Quantile Normalization
def normalize_quantile(df, distribution="normal"):
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
    distribution: str
        Distribution for Quantile normalization. Available distributions are:
        [normal, uniform].

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
        n_quantiles=len(df.T), output_distribution=distribution, 
        random_state=RANDOM_STATE)
    qtransformed_df = pd.DataFrame(
        tf.fit_transform(df.T).T, columns=df.columns, index=df.index)

    # Adjust the scale and location to match original data's mean and sd
    original_mean = df.mean().mean()
    original_std = df.stack().std()
    adj_df = (qtransformed_df * original_std) + original_mean

    # Return the adjusted quantile-transformed data
    return adj_df

## Median Scaling    
def normalize_median_scaling(df, axis:int=0):
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


def normalize_median_scaling_intensity(df, axis=0):
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

def normalize_median(df):
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


def normalize_mean(df):
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

def logfcsq(x, pept, arr, a, b):
    """
    compute squared log fold changes as inner part of delayed normalization
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    pept : int
        index in intensities for given peptide.
    arr : np.array
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
                    x[a]*arr[pept, a, :])/np.nansum(x[b]*arr[pept, b, :])))**2


 
def Hp(x, pept, arr):
    """
    compute H_peptide(x) as middle part of delayed normalization
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    pept : int
        index in intensities for given peptide.
    arr : np.array
        matrix of intensities 
        (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    """
    warnings.warn(
        "This is preliminary code and hasn't been sufficiently tested.", 
        UserWarning)
    
    return np.nansum(
            list(
                map(
                    lambda xx: 0 if arr[
                        pept, xx[0], :].sum() == 0 or arr[pept, xx[1], :].sum() == 0 else logfcsq(
                        x, pept, arr, xx[0], xx[1]), 
                        combinations(range(len(x)), 2))))


def H(x, arr):
    """
    ### NEEDS TO BE TESTED, OL
    ### ADDED PEPT AS INPUT, CD
    ### REMOVED PEPT AS IT IS A RUN VARIABLE, OL
    compute H(x) for all peptides
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    arr : np.array
        matrix of intensities 
        (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    """
    warnings.warn(
        "This is preliminary code and hasn't been sufficiently tested.", 
        UserWarning)
    return [Hp(
            x, pept, arr) for pept in range(
                arr.shape[0])]

def Hsum(x, arr):
    """
    ### NEEDS TO BE TESTED, OL
    ### ADDED PEPT AS INPUT, CD
    ### REMOVED PEPT AS IT IS A RUN VARIABLE, OL
    calculate sum of H functions for every peptide.
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    arr : np.array
        matrix of intensities 
        (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    """
    return np.sum(
        [Hp(x, pept, arr) for pept in range(arr.shape[0])])


def normalize_delayed(arr, method="trf", bounds=(0,np.inf)):
    """
    calculate delayed normalization coefficients for given matrix of intensities 
    as introduced by Cox et al., DOI: 10.1074/mcp.M113.031591 
    
    author: BZ
    
    Parameters
    ----------
    df : np.array
        matrix of intensities 
        (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    method : str, optional
        Algorithm for least squares minimization 
        (see scipy.optimize.least_squares), 
        one of {‘trf’, ‘dogbox’, ‘lm’}. The default is "trf".
    bounds : tuple, optional
        Upper and lower bound for minimization algorithm. 
        The default is (0,np.inf).

    Returns
    -------
    res.x : np.array
        normalization factors for all samples and fractions (if present).

    """
    warnings.warn(
        "This is preliminary code and hasn't been sufficiently tested.", 
        UserWarning)

    #global intensities
    x0 = np.ones(arr.shape[1:]).flatten()
    if method in ["trf","dogbox","lm"]:
        res = least_squares(fun=H, x0=x0, method=method, bounds=bounds, xtol=0.01, args=[arr])  
        # Trust-Region-Reflective Algorithm
    else:
        #bounds=[(0.1,1)]*len(x0)
        bounds=[bounds]*len(x0)
        res=minimize(
            fun=Hsum, x0=x0, method=method, bounds=bounds, 
            options={"disp":True, "maxiter":100*arr.shape[0]}, args=arr)
    res = least_squares(
        fun=H, x0=x0, method=method, bounds=(0, np.inf), xtol=0.01, args=[arr])  
    # Trust-Region-Reflective Algorithm
    res = least_squares(
        fun=H,x0=x0,method="lm",max_nfev=2,ftol=10000,args=[arr]) 
    #Levenberg-Marquardt
    return(res.x)

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

#### Do we want to have a logger return???
    
if __name__ == "__main__":
    
    # Get data from wrapper
    df_fp = sys.argv[0]
    param_fp = sys.argv[1]
    out_dir_fp = sys.argv[2]
    
    # To be modified
    outname = "data_imputed"
    
    # Load parameter_file
    with open('Test_folder/4_IMPUTATION/Module_Parameters.json') as f:
        parameters_dict = json.load(f)
    # Load data
    if parameters_dict["normalization_method"] not in ["delayed"]:
        data_df = pd.read_csv(df_fp,sep="\t",index_col=0)
    else:
        data_arr = np.load(df_fp)
        
    # Load mask and log file
    """
    mask = load_mask() # load mask (in parent directory)
    log = load_log() # load log-file
    
    
    # Helper DataFrames
    log = helper.logger(log, parameters_dict)
    """
    
    # Process data
    if parameters_dict["normalization_method"] == "quantile":
        data_processed_df = normalize_quantile(data_df, parameters_dict["distribution"])
        
    elif parameters_dict["normalization_method"] == "median_scaling":
        data_processed_df = normalize_median_scaling(data_df, parameters_dict["axis"])
        
    elif parameters_dict["normalization_method"] == "median_scaling_intensity":
        data_processed_df = normalize_median_scaling_intensity(data_df, parameters_dict["axis"])
        
    elif parameters_dict["normalization_method"] == "median":
        data_processed_df = normalize_median(data_df)
        
    elif parameters_dict["normalization_method"] == "mean":
        data_processed_df = normalize_mean(data_df)
        
    elif parameters_dict["normalization_method"] == "delayed":
        data_processed_arr = normalize_delayed(data_arr, parameters_dict["method"], parameters_dict["bounds"])

    # Comparison for mask and logger what changed between both files
    """
    # Helper DataFrames
    mask = helper.update_mask(mask)
    log = helper.logger(logf)
    """
    
    # Save data
    if parameters_dict["normalization_method"] not in ["delayed"]:
        data_processed_df.to_csv(os.path.join(out_dir_fp,outname+".tsv"),sep="\t")
    else:
        np.save(os.path.join(out_dir_fp,out_dir_fp+".npy"),data_processed_arr)
    
    