"""
Input: wide-format data
Output: Summarised wide-format data

    - weighted_harmonic_mean: calculates weighted harmonic mean
    - geometric_mean: calculates the geometric mean of an array x
    - topN_summarisation: summarizes peptide intensities to protein intensities with modified topN using pythagorean mean or aligned median on representative precursors for a protein
    - _pairwise_complete_corr: helper for topN_summarisation, calculates Pearson correlation between all rows of X over the columns in which both rows are observed
    - aligned_pooling: aligned median for topNsummarisation, pools log2 intensities of precursors of one protein to a protein profile after aligning the precursors
    - (transform_profiles_to_longformat_table: reshapes wide format profiles with IDs as rows, samples as columns and abundances as values to long-format with ID, sample and values as columns)
    - precursor_sum: sums precursor intensities to obtain protein quantities

"""

# TODO
# Adjust doc-strings

#general
import sys
import pandas as pd
import numpy as np
import time
#from pathlib import Path
import os
import warnings

#own scripts
import module_helper as helper 

#import warnings
#warnings.filterwarnings("ignore")

def weighted_harmonic_mean(x,w=None):
    """
    calculate weighted harmonic mean.
    
    author: BZ
    
    Parameters
    ----------
    x : np.array or list
        Array of floats for which the weighted harmonic mean should be calculated.
    w : np.array or list, optional
        weights associated with values in x for the calculation of the weighted harmonic mean. The default is None, resulting in no weightingfor x.

    Returns
    -------
    float : weighted harmonic mean of x

    """
    x=np.array(x)
    if w is not None:
        w=np.array(w)
    else:
        w=np.ones(len(x))
    return(np.sum(w)/np.sum(w/x))
    


def geometric_mean(x):
    """
    Calculate the geometric mean of an array x.
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        Array of floats.

    Returns
    -------
    float

    """
    x=np.asarray(x,dtype=float)
    # in log space to avoid overflow of the product; 0 if any value is 0
    if np.any(x<=0):
        return(0.0)
    return(float(np.exp(np.mean(np.log(x)))))
    


def topN_summarization(data_df, id_mapping_df, filter_mask, module_output, top_n_peptides=3,weights=None,pooling="aligned",min_peptides=1, seq_col="Precursor.ID", protein_col = "Protein.Group"):
    """
    summarize peptide intensities to protein intensities
    top_n_peptides most abundant, highly correlated and least imputed precursors over samples are summarized to a protein intensity
    #ranks precursors on number of measured peptides (not imputed ones)
    #ranks precursors on abundance
    #ranks precursors on correlation
    all three ranks are normalized, their harmonic mean is calculated and the average as defined by pooling over top_n_peptides precursors is calculated to retrieve the protein intensity.
    The selection of the top_n_peptides precursors is global, i.e. the same precursors are used in all samples.
    All three scores use measured values only (present in data_df and not marked as imputed in filter_mask), so the selection does not depend on the imputation:
        abundance: sum of the linear intensities of measured values over samples, normalized to the maximum within the protein group.
        correlation: Pearson correlation of measured log2 intensities with every other precursor of the protein group (pairwise complete, at least 3 shared samples),
            correlations below 0.01 set to 0.01, combined by the harmonic mean and normalized to the maximum within the protein group.
        measured values: number of measured values, normalized to the maximum within the protein group.
    pooling is how to summarize selected precursors to protein group profiles:
        "mean", "hmean", "gmean": arithmetic, harmonic or geometric mean of the linear intensities of the selected precursors per sample.
            Precursors of a protein differ in intensity by orders of magnitude, so "mean" is dominated by the most and "hmean" by the least intense selected precursor.
        "aligned": the selected precursors are aligned before pooling (log space): every precursor is shifted by its offset, the median over the other
            selected precursors of the median log-ratio over samples in which both were observed (as in graphbased_imputation). The protein value per sample
            is the median of the aligned precursors observed in that sample plus the mean offset of the precursors used, i.e. the absolute level corresponds to
            the geometric mean of the selected precursors. If all selected precursors are complete and parallel, "aligned" equals "gmean".
            Fold changes between samples do not depend on which precursors are selected (for parallel precursors); the absolute level does, but only
            through the mean of the selected precursors, not dominated by the least ("hmean") or most ("mean") intense one. The median makes the profile
            robust against single interfered precursor values. Missing precursor values are skipped (no imputation required), and a protein value is
            missing only if none of the selected precursors was observed in that sample.

    author: BZ, Claude Opus 5.5

    Parameters
    ----------
    data_df : pd.DataFrame (old: peptide_profiles)
        rows equal precursor ids, columns correspond to samples and values in the dataframe to log2-transformed intensities (missing values as NaN).
        data_df is not changed.
    id_mapping_df: pd.DataFrame
        information on mapping from precursor to protein (group) level. rows denote individual precursors
    filter_mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The integer part encodes detection
        or filtering (0 not detected, 1 detected, 2+ filter), a non-zero fractional part marks an imputed value (optionally encoding the imputation method).
    module_output: str
        log string to be filled.
    top_n_peptides : int, optional
        Number of precursors to summarize to protein intensities. The default is 3.
    weights : iterable, optional
        Weights for calculation of harmonic mean over the scores. The default is None, assigning equal weights.
    pooling : str, optional
        one of {"mean","hmean","gmean","aligned"} to summarize selected precursors to proteins (see above). The default is "aligned".
        "mean", "hmean" and "gmean" pool missing values as 0 and therefore require complete (imputed) data; on data with missing values a ###WARNING### is
        written to module_output. "aligned" handles missing values.
    min_peptides : int, optional
        require at least min_peptides precursors of the protein group in data_df to calculate a protein intensity, else the protein group is not
        included in pg_profiles (precursors marked with Summarisation = -1 in id_mapping_df_c). The default is 1.
    seq_col : str
        name of index level containing precursor Ids. The default is 'Precursor.Id'.
    protein_col : str
        name of index level containing protein Ids. The default is 'Protein.Group'.

    Returns
    -------
    pg_profiles : pd.DataFrame
        rows equal protein groups, columns correspond to samples and values in the dataframe to log2-transformed intensities (missing values as NaN).
    id_mapping_df_c: pd.DataFrame
        copy of id_mapping_df with index set to seq_col and column "Summarisation" (1 precursor used, 0 not used, -1 protein group below min_peptides).
        id_mapping_df itself is not changed.
    module_output: str
        log string to be filled.

    """
    
    if pooling not in ["mean","hmean","gmean","aligned"]:
        raise ValueError("pooling must be one of 'mean', 'hmean', 'gmean', 'aligned'")

    n_values_partially_missing, n_values_pooled = 0, 0 # for the warning on missing values with pooling mean/hmean/gmean

    if weights is None:
        weights=[1,1,1]

    # log2 data for pooling="aligned"; pooling="mean","hmean","gmean" on linear intensities with missing values as 0
    log_df = data_df.astype(float)
    data_df = helper.unlog_transform(log_df).replace(np.nan,0)
    id_mapping_df_c = id_mapping_df.copy()

    module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tSummarizing precursor quantities to Protein.Group intensities using own method...\n"
    module_output += "\t\t\t\t- summarize top "+str(top_n_peptides)+" precursors by their "+pooling+"\n"
    module_output += "\t\t\t\t- for the summarization precursors were selected by high abundance, correlation and low number of imputed values with following weights: "+", ".join([str(x) for x in weights])+"\n"

    # measured values: present in data_df and not imputed (non-zero fractional part in filter_mask), used for the precursor selection
    mask_values = filter_mask.reindex(index=log_df.index, columns=log_df.columns).astype(float).fillna(0).to_numpy()
    imputed = ~np.isclose(np.mod(mask_values,1), 0)
    measured_log_df = log_df.where(log_df.notna().to_numpy() & ~imputed)
    
    # Check if columns are in id_mapping_df_c
    if protein_col not in id_mapping_df_c:
        module_output += "###ERROR### Group column for summariaztion not found"
        return log_df, id_mapping_df_c, module_output
    
    if seq_col not in id_mapping_df_c:
        module_output += "###ERROR### Sequence column for summariaztion not found"
        return log_df, id_mapping_df_c, module_output
    
    # numpy arrays (precursors x samples) for the per protein group calculations
    L_all = log_df.to_numpy(dtype=float)                       # log2 intensities, missing as NaN (pooling="aligned")
    lin_all = data_df.to_numpy(dtype=float)                    # linear intensities, missing as 0 (pooling="mean","hmean","gmean")
    M_all = measured_log_df.to_numpy(dtype=float)              # measured log2 intensities, missing or imputed as NaN (selection)
    Mlin_all = helper.unlog_transform(measured_log_df).to_numpy(dtype=float)
    prec_ids = log_df.index.astype(str).to_numpy()
    
    # rows of data_df per protein group, in the order of id_mapping_df_c, restricted to precursors in data_df
    if not log_df.index.is_unique:
        raise ValueError("Precursor ids (index of data_df) must be unique")
    rows_all = log_df.index.get_indexer(id_mapping_df_c[seq_col])
    in_data = rows_all >= 0
    pg_rows = pd.Series(rows_all[in_data]).groupby(id_mapping_df_c[protein_col].to_numpy()[in_data], sort=False).apply(np.array).to_dict()
    pgs_list = id_mapping_df_c[protein_col].unique() 
    
    id_mapping_df_c.set_index(seq_col,inplace=True)
    w = np.asarray(weights, dtype=float)
    summarisation_flags = {} # precursor id -> 1 (used) or -1 (protein group below min_peptides), applied after the loop
    pg_names, pg_profiles = [], []
    
    #for every protein group
    for pg in pgs_list:
        
        if pg not in pg_rows:
            continue
        r = pg_rows[pg]
        n_precs = len(r)
        
        ##on weighted harmonic mean between normalized abundance, correlation and amount of measured samples, all on measured values only
        M = M_all[r]
        observed = ~np.isnan(M)
        # abundance: sum of the linear intensities of measured values
        sums = np.nansum(Mlin_all[r], axis=1)
        # correlation: pairwise complete Pearson correlation (at least 3 shared samples) with every other precursor, clipped at 0.01, harmonic mean
        corrs = _pairwise_complete_corr(M, observed, min_periods=3)
        np.fill_diagonal(corrs, np.nan)
        has_corr = ~np.isnan(corrs)
        n_corr = has_corr.sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            corrs_hmean = np.where(n_corr > 0, n_corr / np.where(has_corr, 1.0/np.clip(np.where(has_corr, corrs, 1.0), 0.01, None), 0.0).sum(axis=1), 0.01)
        # number of measured values
        measured_precs = observed.sum(axis=1)
        
        scores = np.column_stack([(sums/sums.max()) if sums.max()>0 else np.zeros(n_precs),
                                  corrs_hmean/np.max(corrs_hmean),
                                  (measured_precs/measured_precs.max()) if measured_precs.max()>0 else np.zeros(n_precs)])
        with np.errstate(divide="ignore"):
            hmean = np.sum(w)/np.sum(w/scores, axis=1)
        
        # ranking for the selection: score, then measured abundance, then precursor ID (reproducible regardless of input order)
        order = np.lexsort((prec_ids[r], -sums, -np.round(hmean, 12)))
        selected = r[order[:top_n_peptides]]
        
        if n_precs>=min_peptides:
            for prec in prec_ids[selected]:
                summarisation_flags[prec] = 1
            if pooling!="aligned":
                # protein x sample values in which some, but not all, selected precursors are missing (pooled as 0)
                sel_missing = np.isnan(L_all[selected])
                n_values_partially_missing += int((sel_missing.any(axis=0) & ~sel_missing.all(axis=0)).sum())
                n_values_pooled += sel_missing.shape[1]
            P = lin_all[selected]
            if pooling=="mean":
                profile = P.mean(axis=0)
            elif pooling=="hmean":
                with np.errstate(divide="ignore"):
                    profile = P.shape[0]/np.sum(1.0/P, axis=0)
            elif pooling=="gmean":
                with np.errstate(divide="ignore"):
                    profile = np.where(np.any(P<=0, axis=0), 0.0, np.exp(np.mean(np.log(np.where(P>0, P, 1.0)), axis=0)))
            elif pooling=="aligned":
                profile = aligned_pooling(L_all[selected])
            pg_names.append(pg)
            pg_profiles.append(profile)
        else:
            for prec in prec_ids[selected]:
                summarisation_flags[prec] = -1
    
    # Summarisation column: 1 precursor used, 0 not used, -1 protein group below min_peptides
    id_mapping_df_c["Summarisation"] = pd.Series(id_mapping_df_c.index.astype(str), index=id_mapping_df_c.index).map(summarisation_flags).fillna(0).astype(int)
    pg_profiles_df = pd.DataFrame(np.vstack(pg_profiles) if pg_profiles else np.empty((0, log_df.shape[1])), index=pg_names, columns=log_df.columns)

    if n_values_partially_missing > 0:
        module_output += ("###WARNING### Summarisation with pooling='"+pooling+"' on data with missing values: in "+str(n_values_partially_missing)+
                          " of "+str(n_values_pooled)+" protein x sample values ("+str(round(100*n_values_partially_missing/n_values_pooled,1))+
                          "%) at least one selected precursor is missing and pooled as 0, which biases the protein value"+
                          (" (hmean/gmean: value lost)" if pooling in ["hmean","gmean"] else " (mean: value too low)")+
                          ". Use pooling='aligned' for data with missing values or impute before summarisation.\n")
    pg_profiles_df = pg_profiles_df.astype(float)
    if pooling!="aligned":
        # linear protein intensities to log2, 0 (no intensity) as missing value
        pg_profiles_df = np.log2(pg_profiles_df.where(pg_profiles_df>0))
    pg_profiles_df.index.name="Protein.Group"

    return(pg_profiles_df,id_mapping_df_c,module_output)

def _pairwise_complete_corr(X, observed=None, min_periods=3):
    """
    Pearson correlation between all rows of X over the columns in which both rows are observed (as pd.DataFrame.T.corr(min_periods)),
    computed with matrix products instead of a loop over pairs.

    author: BZ, Claude Opus 5.5

    Parameters
    ----------
    X : np.array
        rows (e.g. precursors) x columns (e.g. samples), missing values as NaN.
    observed : np.array of bool, optional
        ~np.isnan(X), passed if already computed.
    min_periods : int, optional
        minimum number of shared observed columns, else NaN. The default is 3.

    Returns
    -------
    np.array
        rows x rows correlation matrix, NaN if fewer than min_periods shared columns or zero variance in the shared columns.
    """
    if observed is None:
        observed = ~np.isnan(X)
    V = observed.astype(float)
    # center every row by its mean for numerical stability (Pearson correlation is invariant to shifts)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning) #all-NaN rows
        Xc = np.where(observed, X - np.nanmean(X, axis=1, keepdims=True), 0.0)
    n = V @ V.T
    sx, sy = Xc @ V.T, V @ Xc.T
    sxx, syy = (Xc*Xc) @ V.T, V @ (Xc*Xc).T
    sxy = Xc @ Xc.T
    with np.errstate(divide="ignore", invalid="ignore"):
        cov = sxy - sx*sy/n
        var_x = sxx - sx*sx/n
        var_y = syy - sy*sy/n
        denom = np.sqrt(var_x*var_y)
        corr = np.where(denom > 0, cov/denom, np.nan)
    corr[n < min_periods] = np.nan
    return np.clip(corr, -1.0, 1.0)


def aligned_pooling(log_prec_df, max_chunk_elements=5*10**6):
    """
    Pools log2 intensities of precursors of one protein to a protein profile after aligning the precursors (see topN_summarization, pooling="aligned").
    The offset of precursor p is the median over all precursors q (including p itself, ratio 0) of the median of (log p - log q) over samples in which
    both were observed. Including p itself makes the offsets consistent: for log I[p,s] = a[p] + b[s] the offset is a[p] - median(a) for every p.
    The protein value per sample is the median of (log p - offset p) over precursors observed in that sample plus the mean offset of the precursors used.
    Precursors that were never observed together with another precursor have no offset and are not used, unless they are the only precursor.

    author: BZ, Claude Opus 5.5

    Parameters
    ----------
    log_prec_df : pd.DataFrame or np.array
        rows equal precursors of one protein, columns samples, values log2 intensities (missing values as NaN).
    max_chunk_elements : int, optional
        Maximum number of elements of the intermediate array of pairwise log-ratios (precursors x precursors x samples) per chunk of precursors,
        limits memory use for proteins with many precursors and samples. The default is 5*10**6.

    Returns
    -------
    np.array
        log2 protein intensity per sample, NaN if no precursor with an offset was observed in that sample.
    """
    L = np.asarray(log_prec_df, dtype=float)
    if L.shape[0] == 1:
        return L[0]
    n_precs, n_samples = L.shape
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning) #all-NaN slices
        # pairwise median log-ratios over co-observed samples, R[p,q] = median_s(L[p,s] - L[q,s]), in chunks of precursors p
        R = np.empty((n_precs, n_precs))
        chunk_size = max(1, int(max_chunk_elements//(n_precs*n_samples)))
        for start in range(0, n_precs, chunk_size):
            R[start:start+chunk_size] = np.nanmedian(L[start:start+chunk_size, None, :] - L[None, :, :], axis=2)
        # precursors co-observed with at least one other precursor
        co_observed = ~np.isnan(R)
        np.fill_diagonal(co_observed, False)
        usable = co_observed.any(axis=1)
        offsets = np.nanmedian(R, axis=1)
        if not usable.any():
            return np.full(L.shape[1], np.nan)
        # offsets relative to their mean, so that the absolute level is the mean of the precursors used
        offsets = offsets[usable] - np.mean(offsets[usable])
        aligned = L[usable] - offsets[:, None]
        return np.nanmedian(aligned, axis=0)


# def transform_profiles_to_longformat_table(profiles, id_mapping_df, protein_col, seq_col):
#     """
#     Reshapes wide format profiles with IDs as rows, samples as columns and abundances as values to long-format with ID, sample and values as columns.
    
#     author: BZ
    
#     Parameters
#     ----------
#     profiles : pd.DataFrame
#         rows: precursors, columns: samples
#     id_mapping_df: pd.DataFrame
#         information on mapping from precursor to protein (group) level. rows denote individual precursors
#     seq_col : str
#         name of index level containing precursor Ids. The default is 'Precursor.Id'.
#     protein_col : str
#         name of index level containing protein Ids. The default is 'Protein.Group'.

#     Returns
#     -------
#     lf_df : pd.DataFrame
#         columns: [ID, sample, value]
#     id_mapping_df_c: pd.DataFrame
#         information on mapping from precursor to protein (group) level. rows denote individual precursors. Index set to seq_col.

#     """
#     '''
#     ap_df=pd.DataFrame(columns=["precursor","protein_group","filename","prec_intensity"])
#     ap_df["precursor"]=profiles.index.get_level_values("Precursor.Id").tolist()*len(profiles.columns)
#     ap_df["protein_group"]=profiles.index.get_level_values("Protein.Group").tolist()*len(profiles.columns)
#     ap_df["filename"]=np.repeat(profiles.columns.get_level_values("sample").tolist(),len(profiles))
#     ap_df["prec_intensity"]=profiles.melt()["value"]
#     '''
    
#     #profiles.reset_index(inplace=True)
#     lf_df=profiles.melt(ignore_index=False)
#     lf_df.reset_index(inplace=True)
#     lf_df.columns = [seq_col,"sample","value"]
#     id_mapping_df_c = id_mapping_df.copy().set_index(seq_col)
#     ### Add protein_col
#     lf_df[protein_col] = [id_mapping_df_c.loc[lf_df.loc[index,seq_col],protein_col] for index in lf_df.index]
#     return(lf_df, id_mapping_df_c)

def precursor_sum(data_df, id_mapping_df, module_output, seq_col="Precursor.Id", protein_col="Protein.Group"):
    """
    sums precursor intensities to obtain protein quantities.
    The sum is calculated on linear intensities (input and output are log2-transformed); missing values are skipped,
    proteins without any observed precursor in a sample are missing (NaN).

    author: BZ/OL
    
    Parameters
    ----------
    data_df : pandas.DataFrame
        rows: precursors, columns: samples, values: log2-transformed intensities (missing values as NaN). data_df is not changed.
    id_mapping_df: pd.DataFrame
        information on mapping from precursor to protein (group) level. rows denote individual precursors
    module_output: str
        log string to be filled.
    seq_col : str
        name of index level containing precursor Ids. The default is 'Precursor.Id'.
    protein_col : str
        name of index level containing protein Ids. The default is 'Protein.Group'.

    Returns
    -------
    data_prot_df : pd.DataFrame
        rows equal protein groups, columns correspond to samples and values in the dataframe to log2-transformed intensities.
    module_output: str
        log string to be filled.

    """
    id_mapping_df_c = id_mapping_df.copy().set_index(seq_col)
    data_prot_df = helper.unlog_transform(data_df.astype(float))
    data_prot_df.index = [id_mapping_df_c.loc[index,protein_col] for index in data_prot_df.index]

    data_prot_df = data_prot_df.groupby(by="index",level=0).sum()
    # back to log2, 0 (no observed precursor) as missing value
    data_prot_df = np.log2(data_prot_df.where(data_prot_df>0))

    module_output += "\n"+time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tProteins summarized to Precursors as sum of Precursor intensities.\n"

    return(data_prot_df, module_output)

if __name__ == "__main__":

    module = "SUMMARISATION"
    #set parameters
    if len(sys.argv) <= 1:
        #if script called within python
        ## counter
        ## workdir
        ## data
        ## metadata
        ## log file
        skip_parser = True
        data_df_fp, param_fp, out_dir_fp, meta_df_fp, log_file, idmap_fp, counter, outname = helper.skip_parser_init(module)
    else:
        #if script called from terminal
        skip_parser = False
        data_df_fp = None
        param_fp = None
        out_dir_fp = None
        meta_df_fp = None
        log_file = None
        idmap_fp = None
        counter = None
        outname = None
    exit_module = not skip_parser
    
    data_df_fp, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, mask_description_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = module, 
                                                                                                                                                                           get_metadata=False, 
                                                                                                                                                                           get_id_mapping=True, 
                                                                                                                                                                           get_mask=True, 
                                                                                                                                                                           return_data_fp_only=True,
                                                                                                                                                                           skip_parser=skip_parser,
                                                                                                                                                                           data_df_fp=data_df_fp,
                                                                                                                                                                           param_fp=param_fp,
                                                                                                                                                                           out_dir_fp=out_dir_fp,
                                                                                                                                                                           meta_df_fp=meta_df_fp,
                                                                                                                                                                           log_file=log_file,
                                                                                                                                                                           idmap_fp=idmap_fp,
                                                                                                                                                                           counter=counter,
                                                                                                                                                                           outname=outname)
    
    helper.logger(log_file,time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tLoading input file " + data_df_fp + ".\n")
    try:
        data_df = helper.load_data(data_df_fp)
        helper.logger(log_file,time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tFile loaded successfully.\n")
    except:
        data_df = None
    
    if data_df is None:
        helper.logger(log_file,"###ERROR###: Input data file could not be loaded")
        if exit_module:
            sys.exit("###ERROR###: Input data file could not be loaded")
    
    #if parameters_dict["summarization_method"] in ["mean","gmean","hmean"]:
    #    masked_df = helper.load_mask(data_df_fp)
    #    mask_description_df = helper.load_mask_description(data_df_fp)
    
    # Check if loading was successful, else pass error to log and exit   
    helper.logger(log_file,module_output)
    module_output = ""
    #if metadata_df is None:
    #    sys.exit("###ERROR###: Metadata loading failed.")
    #if "###WARNING###" in module_output:
    #    helper.logger(log_file,"###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])
    #    sys.exit("###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])

    # all summarisation methods work on log2 data; save unlogged precursor data as additional output
    helper.save_data(helper.unlog_transform(data_df), out_dir_fp, outname[:-4]+"_precursors_unlog.tsv")
    helper.logger(log_file,time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tUnlogged precursor data saved to " + os.path.join(out_dir_fp, outname[:-4]+"_precursors_unlog.tsv")+".\n")

    try:
        if parameters_dict["summarization_method"] in ["mean","gmean","hmean","aligned"]:
            data_processed_df, id_mapping_df, module_output = topN_summarization(data_df, id_mapping_df, masked_df, module_output, 
                                                                                 top_n_peptides=parameters_dict["top_n_peptides"],
                                                                                 weights=parameters_dict["weights"],
                                                                                 pooling=parameters_dict["summarization_method"],
                                                                                 min_peptides=parameters_dict["min_peptides"], 
                                                                                 seq_col=parameters_dict["seq_col"], 
                                                                                 protein_col = parameters_dict["protein_col"])
        elif parameters_dict["summarization_method"]=="sum":
            data_processed_df, module_output = precursor_sum(data_df, id_mapping_df, module_output, 
                                                             seq_col=parameters_dict["seq_col"], 
                                                             protein_col = parameters_dict["protein_col"])
        else:
            data_processed_df = data_df
            module_output += f"###WARNING###: Unknown summarization method: {parameters_dict['summarization_method']} — skipping.\n"
    except Exception as e:
        e_type, e_obj, e_tb = sys.exc_info()
        script_name = os.path.split(e_tb.tb_frame.f_code.co_filename)[1]
        helper.logger(log_file, "###ERROR### Module " +script_name + " ran into an error at line : " + str(e_tb.tb_lineno) +"\n\t\t\t\t" + str(e) + "\n")
        if exit_module:
            sys.exit("Summarisation exited with an error.\n")
    #except ValueError:
    #    raise(ValueError)

    if parameters_dict["summarization_method"] in ["mean","gmean","hmean","aligned"]:
        id_mapping_df.to_csv(os.path.join(out_dir_fp,"id_mapping.tsv"),sep="\t")
    
    # summarisation methods return log2 data; save unlogged protein data as additional output
    helper.save_data(helper.unlog_transform(data_processed_df), out_dir_fp, outname[:-4]+"_unlog.tsv")
    module_output += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tUnlogged protein data saved to " + os.path.join(out_dir_fp, outname[:-4]+"_unlog.tsv")+".\n"

    #if parameters_dict["summarization_method"] in ["mean","gmean","hmean"]:
    helper.close_module(data_processed_df, masked_df, mask_description_df, log_file, module_output, out_dir_fp,outname, mask=True)
    #else:
    #    helper.close_module(data_processed_df, masked_df, mask_description_df, log_file, module_output, out_dir_fp,outname, mask=False)
