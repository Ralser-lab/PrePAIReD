"""
Input: wide-format data
Output: Normalized wide-format data
"""

# Imports
import warnings
import numpy as np
import pandas as pd
import os
from scipy.optimize import least_squares, minimize
import time
import sys

from itertools import combinations
from sklearn.preprocessing import QuantileTransformer

import module_helper as helper

####### Not testes yet #########


####################################################################
# Normalization functions
####################################################################  

## Quantile Normalization
def normalize_quantile(df, distribution="normal", seed=42):
    """ 
    Apply quantile normalization to a DataFrame, transforming sample intensities
    to follow a Gaussian distribution with the same grand mean and sd
    as the original data.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        A DataFrame with rows representing precursors and columns representing 
        samples.
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
        n_quantiles=len(df), output_distribution=distribution, 
        random_state=seed)
    qtransformed_df = pd.DataFrame(
        tf.fit_transform(df), columns=df.columns, index=df.index)

    # Adjust the scale and location to match original data's mean and sd
    if distribution == "normal":
        original_mean = df.mean().mean()
        original_std = df.stack().std()
        adj_df = (qtransformed_df * original_std) + original_mean
    elif distribution == "uniform":
        original_min = df.min().min()
        original_max = df.max().max()
        adj_df = qtransformed_df * (original_max-original_min) + original_min

    # Return the adjusted quantile-transformed data
    return adj_df

def normalize_quantile_ranked(df):
    """ 
    Apply quantile normalization (rank ordered quantile normalization).
    Very similiar to uniform distribution.

    Parameters:
    -----------
    df : pandas.DataFrame
        A DataFrame with rows representing precursors and columns representing 
        samples.

    Returns:
    --------
    adj_df : pandas.DataFrame
        A DataFrame with rows representing precursors and columns representing 
        samples.

    Author: SB
    """

    # Rank the values 
    rank_df = df.rank(method="average")
    
    # Compute the average ranked value across all samples
    sorted_mean_df = df.stack().groupby(rank_df.stack().astype(int)).mean()
    
    # Replace the ranked values with the average ranked values
    adj_df = rank_df.stack().astype(int).map(sorted_mean_df).unstack()
    
    # Return the normalized dataframe
    return adj_df

## Median Scaling    
def normalize_median_scaling(df, out_dir_fp, module_output, normalization_on_intersecting_precs=True, nan_intersect_thresh=0.8, exclude_outlier_nfs=False, rel_threshold=2, on_mad=True, min_nf=0.5, max_nf=2, axis:int=0):
    """

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame with rows representing precursors and columns representing 
        samples.
    out_dir_fp : str
        File path of the output directory to save the normalization factors.
    module_output: str
        Log file to be filled.
    normalization_on_intersecting_precs : bool, optional
        if True, performs normalization only on precursors present in nan_intersect_thresh samples.
    nan_intersect_thresh : float, optional
        use for normalization all precursors present in at least nan_intersect_thresh samples. The default is 0.8.
    exclude_outlier_nfs : bool, optional
        switch to True to exclude samples with outlier normalization factors as defined by rel_threshold, on_mad, min_nf, max_nf. The default is False.
    rel_threshold : float, optional
        defines the threshold of acceptable normalization factors by median(normalization factors) +/- rel_threshold times the MAD. The default is 2.
    on_mad : bool, optional
        switch True if filter by median +/- rel_threshold * MAD. The default is True.
    min_nf : float, optional
        if on_mad == False, min_nf defines the minimal accepted normalization factor. The default is 0.5.
    max_nf : float, optional
        if on_mad == False, max_nf defines the maximal accepted normalization factor. The default is 2.
    axis: int, optional
        Axis for calculation of the normalization factors. The default is 0.

    Returns
    -------
    df_norm : pandas.DataFrame
        A DataFrame containing the median scaled data
    module_output, str
        Filled log file.

    """ 
    if normalization_on_intersecting_precs:
        #extract precursors present over all samples
        abs_nan_intersect_thresh=int(nan_intersect_thresh*df.shape[1])
        intersect_precursors=df.replace(0,np.nan).dropna(thresh=abs_nan_intersect_thresh).index
        module_output += "\n"+time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+str(len(intersect_precursors)) + " precursors that are present in " + str(nan_intersect_thresh*100) + "% of all samples are used for normalization across samples"+"\n"
    else:
        intersect_precursors=df.index
    #calculate normalization factors
    normalization_factors=get_normalization_factors(df.loc[intersect_precursors],axis=axis)
    normalization_factors_df=pd.DataFrame(index=df.columns,columns=["NormalizationFactor"],data=normalization_factors)
    normalization_factors_df.to_csv(os.path.join(out_dir_fp,"Normalization_factors.tsv"),sep="\t")
    if exclude_outlier_nfs: #if normalization factor is np.inf, not sufficient precursors where measured in given sample
        df_trim=df
        while True: #as long as normalization factors are not in given range of min_nf and max_nf, exclude outliers and recalculate them
            #exclude samples with outlier normalization factors
            df_trim, module_output=exclude_outlier_normalization_factors(profiles_group=df_trim,normalization_factors=normalization_factors,module_output=module_output,rel_threshold=rel_threshold,on_mad=on_mad,min_nf=min_nf,max_nf=max_nf)
            if normalization_on_intersecting_precs:
                #extract precursors present over all samples
                abs_nan_intersect_thresh=int(nan_intersect_thresh*df_trim.shape[1])
                intersect_precursors=df_trim.replace(0,np.nan).dropna(thresh=abs_nan_intersect_thresh).index
                module_output += "\n"+time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+str(len(intersect_precursors)) + " precursors that are present in " + str(nan_intersect_thresh*100) + "% of all samples are used for normalization across samples"+"\n"
            else:
                intersect_precursors=df_trim.index            #repeat normalization procedure on selected samples
            normalization_factors=get_normalization_factors(df_trim.loc[intersect_precursors],axis=axis)
            if (np.min(normalization_factors)>min_nf) & (np.max(normalization_factors)<max_nf):
                break
        
        # save final
        normalization_factors_df=pd.DataFrame(index=df_trim.columns,columns=["NormalizationFactor"],data=normalization_factors)
        normalization_factors_df.to_csv(os.path.join(out_dir_fp,"Normalization_factors_trim.tsv"),sep="\t")
        df_norm=df_trim*normalization_factors
    else:
        df_norm=df*normalization_factors
        
    return df_norm, module_output

def exclude_outlier_normalization_factors(profiles_group,normalization_factors,module_output,group,rel_threshold=2,on_mad=True,min_nf=0.5,max_nf=2):
    """
    determine samples with outlier normalization factors by median +- rel_threshold * MAD or by specified normalization factors.    
    Parameters
    ----------
    profiles_group : pd.DataFrame
        rows equal protein groups, columns correspond to samples and values in the dataframe to intensities.
    normalization_factors : pd.Series
        Series of normalization factors for profiles_group.columns.
    module_output: str
        Log file to be filled.
    rel_threshold : float, optional
        defines the threshold of acceptable normalization factors by median(normalization factors) +/- rel_threshold times the MAD. The default is 2.
    on_mad : bool, optional
        switch True if filter by median +/- rel_threshold * MAD. The default is True.
    min_nf : float, optional
        if on_mad == False, min_nf defines the minimal accepted normalization factor. The default is 0.5.
    max_nf : float, optional
        if on_mad == False, max_nf defines the maximal accepted normalization factor. The default is 2.

    Returns
    -------
    profiles_group : pd.DataFrame
        rows equal protein groups, multi columns correspond to samples and values in the dataframe to intensities, samples with normalization factors out of range excluded.
    module_output, str
        Filled log file.
        
    """
        
    normalization_factors_df=pd.DataFrame(index=profiles_group.columns,columns=["NormalizationFactor"],data=normalization_factors)
    
    if on_mad:
        max_nf=normalization_factors_df["NormalizationFactor"].median() + rel_threshold * helper.mad(normalization_factors_df["NormalizationFactor"])
        min_nf=normalization_factors_df["NormalizationFactor"].median() - rel_threshold * helper.mad(normalization_factors_df["NormalizationFactor"])
    normalization_factors_df["keep"]=normalization_factors_df.applymap(lambda x: min_nf < x < max_nf)
    #print("Following samples were excluded due to extreme normalization factors in group "+group+":")
    #print(normalization_factors[~normalization_factors["keep"]].index.tolist())
    module_output += "\t\t\t\tFollowing samples were excluded due to extreme normalization factors:\n\t\t\t\t"
    module_output += ", ".join(normalization_factors_df[~normalization_factors_df["keep"]].index.tolist())+"\n\n"
    return [profiles_group.loc[normalization_factors_df.index[normalization_factors_df["keep"]]], module_output]

def get_normalization_factors(df, axis:int=0):
    """
    Perform median scaling normalization.
    
    author: BZ
    
    Parameters
    ----------
    intensities : pandas.DataFrame
        A DataFrame with rows representing precursors and columns representing 
        samples.

    Returns
    -------
    np.array
        Normalization factors for all samples.
        
    """
    # Throw a ValueError if axis is not 1 or 0
    if axis not in [0, 1]:
        raise ValueError("Axis must be either 0 (rows) or 1 (columns)")

    return np.nanmedian(np.nanmedian(df, axis=axis))/np.nanmedian(df,axis=axis)

def normalize_from_peptide_list(data_df,normpeptides_list,id_mapping_df,output_dir_fp,normlabel=None,precursor_col="Precursor.Id",peptide_col="Stripped.Sequence",normalization_on_intersecting_precs=True,nan_intersect_thresh=0.9,exclude_outlier_nfs=False,rel_threshold=2,on_mad=False,min_nf=1/2,max_nf=2,axis=0):
    """
    Perform median scaling from a subset of peptides, e.g. an internal peptide standard.
    
    author: BZ

    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal protein groups, columns correspond to samples and values in the dataframe to intensities.
    normpeptides_list : list
        list of peptides selecting subset of peptides used for normalization.
    id_mapping_df : pd.DataFrame
        DataFrame containing a column named precursor_col matching index in data_df and a column peptide_col matching peptides in normpeptides_list; index is irrelevant.
    out_dir_fp : str
        File path of the output directory to save the normalization factors.
    normlabel : str or None, optional
        a substring selecting only precursors containing normlabel and being part of normpeptides_list. E.g. "-H" would select only heavy labelled precursors from normpeptides_list. The default is None.
    precursor_col : str, optional
        Column name in id_mapping_df specifying column with precursor ids. The default is "Precursor.Id".
    peptide_col : str, optional
        Column name in id_mapping_df specifying column with peptide sequences. The default is "Stripped.Sequence".
    normalization_on_intersecting_precs : bool, optional
        if True, performs normalization only on precursors present in nan_intersect_thresh samples.
    nan_intersect_thresh : float, optional
        use for normalization all precursors present in at least nan_intersect_thresh samples. The default is 0.8.
    exclude_outlier_nfs : bool, optional
        switch to True to exclude samples with outlier normalization factors as defined by rel_threshold, on_mad, min_nf, max_nf. The default is False.
    rel_threshold : float, optional
        defines the threshold of acceptable normalization factors by median(normalization factors) +/- rel_threshold times the MAD. The default is 2.
    on_mad : bool, optional
        switch True if filter by median +/- rel_threshold * MAD. The default is True.
    min_nf : float, optional
        if on_mad == False, min_nf defines the minimal accepted normalization factor. The default is 0.5.
    max_nf : float, optional
        if on_mad == False, max_nf defines the maximal accepted normalization factor. The default is 2.
    axis: int, optional
        Axis for calculation of the normalization factors. The default is 0.

    Returns
    -------
    data_norm_df : pd.DataFrame.
        dataframe containing normalized intensities
    log : str.
        string containing report for logfile.

    """
    #initialize log
    log=""
    # get dataframe with precursors originating from normpeptides_list
    if normlabel is not None:
        normprecursors_df = data_df.loc[list(filter(lambda x: normlabel in x, id_mapping_df[id_mapping_df["Stripped.Sequence"].isin(normpeptides_list)][precursor_col]))]
    else:
        normprecursors_df = data_df.loc[id_mapping_df[id_mapping_df["Stripped.Sequence"].isin(normpeptides_list)][precursor_col]]
    # select a subset of IDs which is present in at least nan_intersect_thresh *100% samples
    if normalization_on_intersecting_precs:
       #extract precursors present over all samples
       abs_nan_intersect_thresh = int(nan_intersect_thresh*normprecursors_df.shape[1])
       intersect_precursors = normprecursors_df.replace(0,np.nan).dropna(thresh=abs_nan_intersect_thresh).index
       log += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+str(len(intersect_precursors)) + " precursors that are present in " + str(nan_intersect_thresh*100) + "% of all samples are used for normalization across samples.\n"
    else:
        intersect_precursors = data_df.index
    nfs = get_normalization_factors(normprecursors_df.loc[intersect_precursors])
    normalization_factors_df=pd.DataFrame(index=data_df.columns,columns=["NormalizationFactor"],data=nfs)
    normalization_factors_df.to_csv(os.path.join(out_dir_fp,"Normalization_factors.tsv"),sep="\t")

    if exclude_outlier_nfs: #if normalization factor is np.inf, not sufficient precursors where measured in given sample
       normprecursors_trim_df = normprecursors_df
       while True: #as long as normalization factors are not in given range of min_nf and max_nf, exclude outliers and recalculate them
           #exclude samples with outlier normalization factors
           normprecursors_trim_df, log = exclude_outlier_normalization_factors(profiles_group=normprecursors_trim_df,normalization_factors=nfs,module_output=log,rel_threshold=rel_threshold,on_mad=on_mad,min_nf=min_nf,max_nf=max_nf)
           if normalization_on_intersecting_precs:
              #extract precursors present over all samples
              abs_nan_intersect_thresh = int(nan_intersect_thresh*normprecursors_trim_df.shape[1])
              intersect_precursors = normprecursors_trim_df.replace(0,np.nan).dropna(thresh=abs_nan_intersect_thresh).index
              log += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+str(len(intersect_precursors)) + " precursors that are present in " + str(nan_intersect_thresh*100) + "% of all samples are used for normalization across samples.\n"
           else:
               intersect_precursors = data_df.index
           #repeat normalization procedure on selected samples
           nfs=get_normalization_factors(normprecursors_trim_df.loc[intersect_precursors],axis=axis)
           if (np.min(nfs)>min_nf) & (np.max(nfs)<max_nf):
               break
       #save final normalization factors
       normalization_factors_df = pd.DataFrame(index=normprecursors_trim_df.columns,columns=["NormalizationFactor"],data=nfs)
       normalization_factors_df.to_csv(os.path.join(out_dir_fp,"Normalization_factors_trim.tsv"),sep="\t")
       data_norm_df = data_df[normprecursors_trim_df.coluns] * nfs
    else:
       data_norm_df = data_df * nfs
       
    return([data_norm_df,log])


def normalize_median(df):
    """
    median normalization 
    substracting median of every column, resulting in column medians of 0.
    
    author: BZ
    
    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame with rows representing precursors and columns representing 
        samples.

    Returns
    -------
    intensities : pandas.DataFrame
        pandas.DataFrame of normalized intensities with column medians of 0 
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
    df : pandas.DataFrame
        A DataFrame with rows representing precursors and columns representing 
        samples.

    Returns
    -------
    intensities : pandas.DataFrame
        pandas.DataFrame of normalized intensities with column means of 0 
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


def get_normalization_factors_delayed(arr, module_output, method="trf", bounds=(0,np.inf)):
    """
    calculate delayed normalization coefficients for given matrix of intensities 
    as introduced by Cox et al., DOI: 10.1074/mcp.M113.031591 
    
    author: BZ
    
    Parameters
    ----------
    arr : np.array
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
    module_output += "### WARNING: This is preliminary code and hasn't been sufficiently tested.\n"

    #global intensities
    x0 = np.ones(arr.shape[1:]).flatten()
    if method in ["trf","dogbox","lm"]:
        res = least_squares(fun=H, x0=x0, method=method, bounds=bounds, xtol=0.01, args=[arr])  
        # Trust-Region-Reflective Algorithm
    else:
        module_output += "Method not found. Did not apply least square minimization.\n"
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
    return(res.x), module_output


def normalize_delayed(arr, module_output, method="trf", bounds=(0,np.inf)):
    """
    Parameters
    ----------
    arr : np.array
        matrix of intensities 
        (1st dim peptides, 2nd dim samples, 3rd dim fractions).
    module_output: str
        Log file to be filled.
    method : str, optional
        Algorithm for least squares minimization 
        (see scipy.optimize.least_squares), 
        one of {‘trf’, ‘dogbox’, ‘lm’}. The default is "trf".
    bounds : tuple, optional
        Upper and lower bound for minimization algorithm. 
        The default is (0,np.inf).

    Returns
    -------
    norm_arr : np.array
        Normalized array.
    module_output, str
        Filled log file.
    """
    
    x0= np.ones(arr.shape[1:]).flatten()
    normalization_factors, module_output = get_normalization_factors_delayed(arr,module_output,bounds=bounds,method=method)
    if (normalization_factors==x0).sum()==len(normalization_factors): #no result found
        module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tNormalizing using delayed normalization.\n"
        module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tMinimization to determine normalization factors ended at initial values\n"
    norm_arr = arr*normalization_factors.T
    return norm_arr, module_output


if __name__ == "__main__":

    if parameters_dict["normalization_method"] in ["delayed"]:
        data_df, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, mask_description_df, out_dir_fp, outname, log_file, counter, data_df_fp = helper.init_module(module = "Normalisation", metadata=True, id_mapping=True, mask=True, return_data_fp=True)
    else:
        data_df, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, mask_description_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = "Normalisation", metadata=True, id_mapping=True, mask=True, return_data_fp=False)
    
    # Load data
    if parameters_dict["normalization_method"] in ["delayed"]:
        data_arr = np.load(data_df_fp)

    try:
        # Process data
        if parameters_dict["normalization_method"] == "quantile":
            module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Performing quantile normalization.\n"
            data_processed_df = normalize_quantile(data_df, parameters_dict["distribution"], seed)
            
        elif parameters_dict["normalization_method"] == "quantile_ranked":
            module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Performing quantile ranked normalization.\n"
            data_processed_df = normalize_quantile_ranked(data_df)
        elif parameters_dict["normalization_method"] == "median_scaling":
            if parameters_dict["on_peptide_list"]:
               # check if file at peptide_list_fp exists, if not exit with an error
               if not os.path.exists(parameters_dict["peptide_list_fp"]):
                   module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"###ERROR### File path to file containing peptides for normalisation is not existent.\n"
                   helper.logger(log_file, module_output)
                   sys.exit("Error: file path to file containing peptides for normalisation is not existent.")
               # read file containing peptide list
               peptide_df = pd.read_csv(parameters_dict["peptide_list_fp"],sep="\t")
               # check if column with peptides and precursors are present in peptide_df
               if parameters_dict["peptide_list_file_peptide_column"] not in peptide_df.columns or parameters_dict["peptide_list_file_precursor_column"] not in peptide_df.columns:
                   module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"###ERROR### File containing peptides for normalisation ("+ parameters_dict["peptide_list_fp"] +") does not contain column specified by peptide_list_file_peptide_column or peptide_list_file_precursor_column.\n"
                   helper.logger(log_file, module_output)
                   sys.exit("Error: File with peptides for normalization does not contain the given column(s).")
               # extract peptide list
               peptide_list = peptide_df[parameters_dict["peptide_list_file_peptide_column"]].tolist()
               module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Performing median scaling normalization from peptides specified in " + parameters_dict["peptide_list_fp"] + ".\n"
               # perform normalisation
               data_processed_df, log = normalize_from_peptide_list(data_df, peptide_list, id_mapping_df, out_dir_fp,
                                                                      normlabel=parameters_dict["normalisation_precursors_label"],
                                                                      precursor_col=parameters_dict["peptide_list_file_precursor_column"],
                                                                      peptide_col=parameters_dict["peptide_list_file_peptide_column"],
                                                                      normalization_on_intersecting_precs=parameters_dict["normalization_on_intersecting_precs"],
                                                                      nan_intersect_thresh=parameters_dict["nan_intersect_thresh"],
                                                                      exclude_outlier_nfs=parameters_dict["exclude_outlier_nfs"],
                                                                      rel_threshold=parameters_dict["rel_threshold"],
                                                                      on_mad=parameters_dict["on_mad"],
                                                                      min_nf=parameters_dict["min_nf"],
                                                                      max_nf=parameters_dict["max_nf"],
                                                                      axis=parameters_dict["axis_median_scaling"])
               module_output += log

            else:
               module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Performing median scaling normalization.\n"
               data_processed_df, module_output = normalize_median_scaling(data_df, out_dir_fp, module_output,
                                                                      normalization_on_intersecting_precs=parameters_dict["normalization_on_intersecting_precs"],
                                                                      nan_intersect_thresh=parameters_dict["nan_intersect_thresh"],
                                                                      exclude_outlier_nfs=parameters_dict["exclude_outlier_nfs"],
                                                                      rel_threshold=parameters_dict["rel_threshold"],
                                                                      on_mad=parameters_dict["on_mad"],
                                                                      min_nf=parameters_dict["min_nf"],
                                                                      max_nf=parameters_dict["max_nf"],
                                                                      axis=parameters_dict["axis_median_scaling"])
    
        elif parameters_dict["normalization_method"] == "median":
            module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Performing median normalization.\n"
            data_processed_df = normalize_median(data_df)
            
        elif parameters_dict["normalization_method"] == "mean":
            module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Performing mean normalization.\n"
            data_processed_df = normalize_mean(data_df)
    
        elif parameters_dict["normalization_method"] == "delayed":
            #### Not properly tested
            module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Performing delayed normalization.\n"
            data_processed_arr, module_output = normalize_delayed(data_arr, module_output, parameters_dict["method"], parameters_dict["bounds"])
        else:
            data_processed_df = data_df
            module_output += "###WARNING###: Normalzation method not implemented. Skipped module."
        
        module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Normalization finished\.n"
    except:
        helper.logger(log_file, "###ERROR###: Module ran into an error.")
        sys.exit(1)

    # Comparison for mask and logger what changed between both files
    """
    # Helper DataFrames
    mask = helper.update_mask(mask)
    """
    
    helper.logger(log_file, module_output=module_output)
    
    # Save data
    if parameters_dict["normalization_method"] not in ["delayed"]:
        helper.close_module(data_processed_df, masked_df, mask_description_df, log_file, module_output, out_dir_fp,outname, mask=True)
    else:
        np.save(os.path.join(out_dir_fp,out_dir_fp.split("/")[-1]+".npy"),data_processed_arr)
        helper.logger(log_file, module_output=module_output)
        # Save mask
        helper.save_mask(masked_df, out_dir_fp, outname)
        helper.save_mask_description(mask_description_df, out_dir_fp, outname)
        helper.logger(log_file, "Module finished as expected.")
        sys.exit(0)
    
    
    
    
