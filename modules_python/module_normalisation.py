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

#statistics and modelling
from scipy.optimize import least_squares, minimize
from itertools import combinations
from sklearn.preprocessing import QuantileTransformer


def log_transform(df, base = 2):
    """ Log scales a dataframe (no missing values)
    Author: SB
    :param df: Intensity dataframe (rows [samples], columns [protein feature])
    :param base: base for log transformation, default is 2
    :return: Log transformed intensity dataframe
    """

    # Check for missing values in matrix
    if df.isna().any().any():
        raise ValueError("Error: DataFrame contains missing values (NaNs). Please handle them before log transformation.")
    
    # Apply log transformation and recover true Zero
    l_df = np.log(df) / np.log(base)
    l_df.replace(-np.inf, 0, inplace=True)

    # Return log transformed df
    return l_df


def quantile_transform(df):
    """ For scaling experiments with diverse signal distribution (due to severe defects in LC-MS acquisition)
    Brute forces samples (rows) to follow gaussian distribution with same parameters (grand mean and std of original data)
    Author: SB
    :param df: Intensity dataframe (rows [samples], columns [protein feature])
    :return: Quantile transformed intensity dataframe
    """

    # Apply quantile normalization on rows (samples)
    transformer = QuantileTransformer(n_quantiles=len(df.T), output_distribution='normal', random_state=0)
    data_qtransformed = pd.DataFrame(transformer.fit_transform(df.T).T, columns=df.columns, index=df.index)

    # Adjust the scale and location to match original data's mean and standard deviation
    original_mean = df.mean().mean()
    original_std = df.stack().std()
    data_adjusted = (data_qtransformed * original_std) + original_mean

    # Return the adjusted quantile-transformed data
    return data_adjusted

def quantile_normalize(df):
    """ Apply biological quantile normalization (rank ordered quantile normalization for Omics) on a dataframe
    Author: SB
    :param df: Intensity dataframe (rows [samples], columns [protein feature])
    :return: Biological quantile normalized intensity dataframe
    """

    # Rank the values 
    df_rank = df.rank(method="average")
    
    # Compute the average ranked value across all samples
    df_sorted_mean = df.stack().groupby(df_rank.stack().astype(int)).mean()
    
    # Replace the ranked values with the average ranked values
    df_normalized = df_rank.stack().astype(int).map(df_sorted_mean).unstack()
    
    # Return the normalized dataframe
    return df_normalized


def z_scale(df): 
    """ Z-scales each feature (column) in a dataframe, i.e. for covariation / comparing across features
    Author: SB
    :param df: Intensity dataframe (rows [samples], columns [protein feature])
    :return: Z-scaled intensity dataframe
    """
    
    #calculate feature distribution parameters (mean and standard deviation)
    dfmean = df.mean(skipna = True)
    dfstd = df.std(skipna = True)

    #center and scale the features
    df_z = (df - dfmean)/dfstd

    #return z-scaled data
    return df_z 

# %%
####################################################################
# log transformation
####################################################################



# applies log transformation to wide format data (before any annotations)
def log_transform_intensities(profiles):
    """
    log2 transform intensities of wide-format precursor or protein intensities over samples.
    Keep in mind that values between 0 and 1 are lower than 0 after log-transformation, values of 0 are set to 0 again after logtransformation (instead of -np.inf) to make use in other programming languages easier.
    
    author: BZ
    
    Parameters
    ----------
    profiles : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities.
    

    Returns
    -------
    profiles : pd.DataFrame
        same shape as input dataframe, but values are log2 transformed.

    """
   
    #replace NaNs with zeros
    profiles.replace(0,np.nan,inplace=True)
    #transform intensities to log2 intensities
    profiles = np.log2(profiles)
    profiles.replace(np.nan,0,inplace=True)
    return(profiles)


# transforms wide format data back from log2-space
def unlog_transform_intensities(profiles):
    """
    transforms intensities of wide-format precursor or protein intensities over samples from log2-space to not-log2-space.
        
    author: BZ
    
    Parameters
    ----------
    profiles : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities.
    
    Returns
    -------
    profiles : pd.DataFrame
        same shape as input dataframe, but values are log2 transformed.

    """

    #replace NaNs with zeros
    profiles.replace(0,np.nan,inplace=True)
    #transform intensities to log2 intensities
    profiles = 2**profiles
    profiles.replace(np.nan,0,inplace=True)
    return(profiles)






####################################################################
# median scaling
####################################################################


    
def median_scaling_normalization(intensities):
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
    return(np.nanmedian(np.nanmedian(intensities,axis=0))/np.nanmedian(intensities,axis=0))



def median_scaling_intensity(intensities,axis=0):
    """
    Obtain median over median intensities (one value) which is used to normalize a test set to the training set.
    
    author: BZ
    
    Parameters
    ----------
    intensities : np.array
        matrix of intensities 
    axis : int
        axis on which to calculate first median, needs to be zero if 1st dimension of intensities is precursors, 2nd samples, one if transposed. Default is 0.
    Returns
    -------
    float
        median over median intensities

    """
    return(np.nanmedian(np.nanmedian(intensities,axis=axis)))


####################################################################
# z-scoring
####################################################################

def median_zscoring(data,axis=0):
    """
    calculate z-scores over samples as (X-median(X))/MAD(X)
    
    author: BZ
    
    Parameters
    ----------
    data : pd.DataFrame
        dataframe with features as columns and samples as rows
    axis : int, optional
        Set to 0 if columns in dataframe are features and samples are rows, set to 1 if rows and columns are transposed. The default is 0.

    Returns
    -------
    zdata : pd.DataFrame
        dataframe with same shape as data and z-scores as values

    """
    mads=data.sub(data.median(axis=axis),axis=1-axis).abs().median(axis=axis)
    zdata=data.sub(data.median(axis=axis),axis=1-axis).div(mads,axis=1-axis)
    #sanity check
    ##(data.iloc[0,0]-np.median(data.iloc[:,0]))/mads[0] - zdata.iloc[0,0] == 0
    ##zdata.median() #approx. 0
    return(zdata)


####################################################################
# functions for delayed normalization (Cox et al., DOI: 10.1074/mcp.M113.031591)
# require more testing, not tested sufficiently for fractionated data, poorly compared to other normalization strategies
####################################################################

def logfcsq(x, peptide, intensities, a, b):
    """
    compute squared log fold changes as inner part of delayed normalization
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    peptide : int
        index in intensities for given peptide.
    intensities : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    a : int
        sample a.
    b : int
        sample b.


    """
    return(abs(np.log2(np.nansum(x[a]*intensities[peptide, a, :])/np.nansum(x[b]*intensities[peptide, b, :])))**2)


# 
def Hp(x, peptide, intensities):
    """
    compute H_peptide(x) as middle part of delayed normalization
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    peptide : int
        index in intensities for given peptide.
    intensities : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples, 3rd dim fractions).


    """
    return(np.nansum(list(map(lambda xx: 0 if intensities[peptide, xx[0], :].sum() == 0 or intensities[peptide, xx[1], :].sum() == 0 else logfcsq(x, peptide, intensities, xx[0], xx[1]), combinations(range(len(x)), 2)))))



def H(x, intensities):
    """
    compute H(x) for all peptides
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    intensities : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples, 3rd dim fractions).


    """
    return([Hp(x, peptide, intensities) for peptide in range(intensities.shape[0])])


def Hsum(x, intensities):
    """
    calculate sum of H functions for every peptide.
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        normalization factors (samples x fractions).
    intensities : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples, 3rd dim fractions).


    """
    return(np.sum([Hp(x, peptide, intensities) for peptide in range(intensities.shape[0])]))


#TODO
#logging inside this function, because more information to be logged than for other normalization strategies. 
#Could be moved outside function, but would change output.
def delayed_normalization(intensities,method="trf",bounds=(0,np.inf),logfilename=None):
    """
    calculate delayed normalization coefficients for given matrix of intensities as introduced by Cox et al., DOI: 10.1074/mcp.M113.031591 
    
    author: BZ
    
    Parameters
    ----------
    intensities : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    method : str, optional
        Algorithm for least squares minimization (see scipy.optimize.least_squares), one of {‘trf’, ‘dogbox’, ‘lm’}. The default is "trf".
    bounds : tuple, optional
        Upper and lower bound for minimization algorithm. The default is (0,np.inf).

    Returns
    -------
    res.x : np.array
        normalization factors for all samples and fractions (if present).

    """
    #global intensities
    x0 = np.ones(intensities.shape[1:]).flatten()
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tBegin delayed normalization...\n")
        logfile.close()
    if method in ["trf","dogbox","lm"]:
        res = least_squares(fun=H, x0=x0, method=method, bounds=bounds, xtol=0.01, args=[intensities])  # Trust-Region-Reflective Algorithm
    else:
        #bounds=[(0.1,1)]*len(x0)
        bounds=[bounds]*len(x0)
        res=minimize(fun=Hsum, x0=x0, method=method, bounds=bounds, options={"disp":True, "maxiter":100*intensities.shape[0]}, args=intensities)
    #res = least_squares(fun=H, x0=x0, method=method, bounds=(0, np.inf), xtol=0.01, args=[intensities])  # Trust-Region-Reflective Algorithm
    # res=least_squares(fun=H,x0=x0,method="lm",max_nfev=2,ftol=10000,args=[intensities]) Levenberg-Marquardt
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tEnd normalization\n")
        logfile.write("\t\t\t\tReason for minimization termination: "+res.message+"\n")
        try:
            logfile.write("\t\t\t\t"+res.nfev + "\n")
        except:
            pass
        logfile.close()
    return(res.x)


####################################################################
# functions used for ratio-based analysis mainly (still normalization functions)
# e.g. normalization of sample fold changes, normalization to a mean
# might be excluded for now
####################################################################



def median_normalization(intensities):
    """
    median normalization (substracting median of every column, resulting in column medians of 0).
    
    author: BZ
    
    Parameters
    ----------
    intensities : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples)

    Returns
    -------
    intensities : np.array
        matrix of normalized intensities with column medians of 0 (1st dim peptides, 2nd dim samples)

    """
    return(intensities.sub(intensities.median(axis=0),axis=1))


def mean_normalization(intensities):
    """
    mean normalization (substracting mean of every column, resulting in column means of 0).
    
    author: BZ
    
    Parameters
    ----------
    intensities : np.array
        matrix of intensities (1st dim peptides, 2nd dim samples)

    Returns
    -------
    intensities : np.array
        matrix of normalized intensities with column means of 0 (1st dim peptides, 2nd dim samples)

    """
    return(intensities.sub(intensities.mean(axis=0),axis=1))


def calculate_mean_over_refs(df,refstring=r'^BY4741.*'):
    """
    Calculate a mean sample over reference samples specified via refstring and exclude original samples specified via refstring.
    
    author: BZ
    
    Parameters
    ----------
    df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.
    refstring : str, optional
        regular expression used to identify reference samples. The default is r'^BY4741.*'.

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


def calculate_log2FC_to_ref(df,ref,level=None):
    """
    Calculates log2 fold changes to a reference sample specified by ref.
    
    author: BZ
    
    Parameters
    ----------
    df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.
    ref : str
        column name of reference sample to which fold changes should be calculated.
    level : str, optional
        For multi-column df please provide a level in which ref is found. The default is None.

    Returns
    -------
    log2fc_df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions, values are log2 fold changes to reference column.

    """
    if level is None:
        log2fc_df=np.log2(df.div(df[ref],axis=0))
    else:
        log2fc_df=np.log2(df.div(df.xs(ref,axis=1,level=level).iloc[:,0].tolist(),axis=0))
    return(log2fc_df)


def median_log2foldchange(intensities):
    """
    calculate for every id log2 fold change to median intensity over all samples
    
    author: BZ
    
    Parameters
    ----------
    intensities : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.

    Returns
    -------
    pd.DataFrame

    """
    intensities=intensities.replace(0,np.nan)
    log2fcs=np.log2(intensities.div(intensities.median(axis=1),axis=0))
    log2fcs=log2fcs.replace(np.nan,0)
    return(log2fcs)
    

def mean_log2foldchange(intensities):
    """
    calculate for every id log2 fold change to mean intensity over all samples
    
    author: BZ
    
    Parameters
    ----------
    intensities : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.

    Returns
    -------
    pd.DataFrame

    """
    intensities=intensities.replace(0,np.nan)
    log2fcs=np.log2(intensities.div(intensities.mean(axis=1),axis=0))
    log2fcs=log2fcs.replace(np.nan,0)
    return(log2fcs)


def mad(x):
    """
    calculate median absolute deviation of a vector x
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        Array with float values.

    Returns
    -------
    float

    """
    return(np.nanmedian(np.abs(x-np.nanmedian(x))))


####################################################################
#logging EXAMPLE in main
####################################################################

def main(intensities,logfilename=None):
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tNormalizing using median scaling normalization.\n")
        logfile.close()
    median_scaling_intensity(intensities,axis=0)
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tFinished normalization.\n")
        logfile.close()
    
=======
    mean = df.mean(skipna = True)
    #calculate variance
    var = df.var(skipna = True)
    #standardize df
    df_z = (df - mean)/np.sqrt(var)
    #return z-normalized data
    return df_z 

def main():
    return example()
