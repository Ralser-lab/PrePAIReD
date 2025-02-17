# %%
"""
Input: wide-format data
Output: Normalized wide-format data
"""

#Import
import numpy as np
import pandas as pd
import re
import time #used for logging
import warnings

#statistics and modelling
from scipy.optimize import least_squares, minimize
from itertools import combinations
from sklearn.preprocessing import QuantileTransformer

# Global 
RANDOM_STATE = 42

# %%
####################################################################
# Helper functions
####################################################################

def log_transform(
        df, base:int=2, 
        add_constant=False, constant=1):
    """
    Applies a logarithmic transformation to a DataFrame with the specified base.

    Authors: SB, BZ, CD

    Parameters:
    -----------
    df : pandas.DataFrame
        The input DataFrame containing numerical values to be transformed.
        
    base : int (default=2)
        The base of the logarithm. 
        
    add_constant : bool, optional (default=False)
        If True, a constant value will be added to all elements of the DataFrame 
        before applying the logarithm transformation. This can be useful if your 
        data contains zero or negative values (since log of zero or negative is 
        undefined).
        
    constant : int, float, optional (default=1)
        The constant value to be added to each element of the DataFrame, used 
        only if `add_constant=True`.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame with the log-transformed values. The transformed values 
        will have the same shape as the input `df`.
        Missing values (NaN) will remain NaN in the output DataFrame.

    Raises:
    -------
    ValueError:
        If the DataFrame contains non-numeric values or if the base is 
        non-positive or equal to 1.
    """

    # Check for invalid base values
    if base <= 0 or base == 1:
        raise ValueError("Base must be greater than 0 and not equal to 1")
    # Get NaNs
    nan_mask = df.isna()
    # Add constant
    if add_constant:
        df = df + constant
    # Log transform
    df = np.log(df) / np.log(base)
    # Replace Nan's
    df[nan_mask] = np.nan
    return df


def unlog_transform(df, base=2):
    """
    Transforms intensities of wide-format precursor or protein intensities 
    over samples from log_{Base}-space to not-log_{Base}-space.
        
    Author: BZ
    
    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame where rows could represent precursor or protein IDs and 
        columns represent samples. 

    base : int, optional (default=2)
        The base of the logarithmic transformation. By default, it is set to 2 
        for log2-space. You can adjust this to match the base of the original log 
        transformation (e.g., 10 for log10).
    
    Returns
    -------
    df: pd.DataFrame
        A DataFrame of the same shape as the input, but with values transformed 
        back to the original scale by applying the exponentiation with the 
        specified base.

    """
    df = base**df
    return(df)

def mad(x):
    """
    Xalculate median absolute deviation of a vector x
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        Array with float values.

    Returns
    -------
    float 

    """
    return np.nanmedian(np.abs(x-np.nanmedian(x)))

"""
BZ: 
functions used for ratio-based analysis mainly (still normalization functions)
e.g. normalization of sample fold changes, normalization to a mean
might be excluded for now
"""

def calculate_mean_over_refs(df,refstring=r'^BY4741.*'):
    """
    Calculate a mean sample over reference samples specified via refstring and 
    exclude original samples specified via refstring.
    
    author: BZ
    
    Parameters
    ----------
    df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.
    refstring : str, optional
        regular expression used to identify reference samples. 
        The default is r'^BY4741.*'.

    Returns
    -------
    df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.

    """
    #define reference columns
    r=re.compile(refstring)
    refs=[m.group(0) for m in map(r.match,df.columns) if m is not None]
    #calculate mean over reference
    df["ref"]=df[refs].mean(axis=1)
    #remove ref strains
    df.drop(refs,axis=1,inplace=True)
    return(df)


def calculate_log2FC_to_ref(df,ref ,level=None):
    """
    Calculates log2 fold changes to a reference sample specified by ref.
    
    author: BZ
    
    Parameters
    ----------
    df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.
    ref : str
        column name of reference sample to which fold changes should be 
        calculated.
    level : str, optional
        For multi-column df please provide a level in which ref is found. 
        The default is None.

    Returns
    -------
    log2fc_df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions, 
        values are log2 fold changes to reference column.

    """
    if level is None:
        log2fc_df = np.log2(df.div(df[ref], axis=0))
    else:
        log2fc_df = np.log2(
            df.div(df.xs(ref, axis=1,level=level).iloc[:,0].tolist(), axis=0))
    return log2fc_df 


def calculate_log2foldchange(df, robust=True):
    """
    Calculate the log2 fold change for each protein/precursor relative to the 
    median (or mean) intensity across all samples.
    
    author: BZ
    
    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame where rows represent proteins/precursors, and columns 
        represent samples/runs/fractions. The values in the DataFrame are 
        intensity measurements.

    robust : bool, optional (default=True)
        If True, the log2 fold change is calculated relative to the median 
        intensity across samples for each protein. If False, the fold change 
        is calculated relative to the mean intensity across samples.

    Returns
    -------
    pd.DataFrame
        A DataFrame with the same shape as the input, where each value is the 
        log2 fold change of the corresponding protein/precursor relative to 
        the median (or mean) intensity across all samples. Zero values are 
        replaced with `NaN` during calculation and restored to 0 in the result.

    """
    df.replace(0, np.nan, inplace=True)
    if robust:
        log2fcs = np.log2(df.div(df.median(axis=1),axis=0))
    else:
        log2fcs = np.log2(df.div(df.mean(axis=1),axis=0))
    
    log2fcs = log2fcs.replace(np.nan,0)
    return log2fcs 
    

####################################################################
# z-scoring
####################################################################

def zscore_transform(df, robust=False, axis:int=0):
    """

    Calculates z-scores based on mean and sd or median and median absolute 
    deviation (MAD) on the preferred axis (0:rows, 1: columns).

    if robust == True use: z-scores over samples as (X-median(X))/MAD(X)
    
    author: BZ, SB, CD
    
    Parameters
    ----------
    data : pd.DataFrame
        dataframe with features as columns and samples as rows
    axis : int, optional
        Set to 0 if columns in dataframe are features and samples are rows, set 
        to 1 if rows and columns are transposed. The default is 0.

    Returns
    -------
    zdata : pd.DataFrame
        dataframe with same shape as data and z-scores as values

    """

    # Throw a ValueError if axis is not 1 or 0
    if axis not in [0, 1]:
        raise ValueError("Axis must be either 0 (rows) or 1 (columns)")

    if robust:
        mads=df.sub(df.median(axis=axis),axis=1-axis).abs().median(axis=axis)
        df_z=df.sub(df.median(axis=axis),axis=1-axis).div(mads,axis=1-axis)
    else:
        means = df.mean(axis=axis)
        sds = df.std(axis=axis)
        df_z = df.sub(means,axis=1-axis).div(sds, axis=1-axis)

    return(df_z)

# %%

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

    return(
        np.nanmedian(
            np.nanmedian(df, axis=axis))/np.nanmedian(df,axis=axis))


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
    return(df.sub(df.median(axis=0),axis=1))


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
    return(df.sub(df.mean(axis=0),axis=1))



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
    return[Hp(
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


#TODO
#logging inside this function, because more information to be logged than for 
# other normalization strategies. 
# Could be moved outside function, but would change output.
def delayed_normalization(
        df, method="trf", bounds=(0,np.inf), logfilename=None):
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
    x0 = np.ones(df.shape[1:]).flatten()
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(
            time.strftime(
                "%Y-%m-%d%t%H:%M:%S")+"\tBegin delayed normalization...\n")
        logfile.close()
    if method in ["trf","dogbox","lm"]:
        res = least_squares(
            fun=H, x0=x0, method=method, bounds=bounds, xtol=0.01, args=[df])  # Trust-Region-Reflective Algorithm
    else:
        #bounds=[(0.1,1)]*len(x0)
        bounds=[bounds]*len(x0)
        res=minimize(
            fun=Hsum, x0=x0, method=method, bounds=bounds, 
            options={"disp":True, "maxiter":100*df.shape[0]}, args=df)
    #res = least_squares(fun=H, x0=x0, method=method, bounds=(0, np.inf), xtol=0.01, args=[intensities])  # Trust-Region-Reflective Algorithm
    # res=least_squares(fun=H,x0=x0,method="lm",max_nfev=2,ftol=10000,args=[intensities]) Levenberg-Marquardt
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(
            time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tEnd normalization\n")
        logfile.write(
            "\t\t\t\tReason for minimization termination: "+res.message+"\n")
        try:
            logfile.write("\t\t\t\t"+res.nfev + "\n")
        except:
            pass
        logfile.close()
    return(res.x)


####################################################################
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