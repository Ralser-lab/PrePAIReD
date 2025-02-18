"""
Input: X
Output: X
"""

import re
import numpy as np
import pandas as pd


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

    return df_z 

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


def unlog_transform(
        df, base=2):

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

def calculate_mean_over_refs(
        df, refstring=r'^BY4741.*'):
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
    return df


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
    
