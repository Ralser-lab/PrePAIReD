"""
Input: wide-format data
Output: Filtered wide-format data

functions:
    - filter_samples_based_on_metadata: Removes samples of poor quality, based on flags in a designated column in annotated metadata dataframe.
    - filter_samples_by_detection_proportion: Removes samples that contain missing values more than a specified or automated threshold value (missing proportion).
    - filter_precursors_by_detection_proportion: Removes precursors that are detected less than a specified or automated threshold value (missing proportion).
    - filter_precursors_by_exclusionlist: Removes precursors that are specified in a user-provided list (e.g. contaminants) or precursors that map to ids in respective list
    - filter_proteins_by_id_count: Removes precursors or peptides mapping to proteins identified by less than a minimum number of unique precursors or peptides. 
    - filter_precursors_by_metadata: selects precursors with aid of id_mapping_df, e.g. only proteotypic precursors
    - filter_precursors_by_peptide_labels: Removes precursors that contain specific substrings (e.g., SILAC labels like "(SILAC-K-H)").
    - replace_labelled_precursors: Finds and replaces specific labeling substrings (e.g., "(SILAC-K-L)") in precursor IDs with another label (e.g., "K"). 
    - condition_on_detection_proportion: Calculates % missing proportion by conditioning on missingness on either sample or precursor dimension.
    - calculate_optimum_threshold: Automatically determines the most appropriate detection threshold value (for precursor or sample by % missing proportion) based on sampling statistics.
    - update_mask: Updates mask_df with a unique identifier after a specific wide filtering step.

#???
exclusion of samples by outlier normalisation factors currently included in module_normalisation, but without tracking in sample_filters_df

"""

import numpy as np
import pandas as pd
import time
import sys
import os
import scipy.stats as stats
import module_helper as helper

pd.set_option('future.no_silent_downcasting', True)

##############################################################################
# Sample Thresholding Functions
##############################################################################

def filter_samples_based_on_metadata(data_df, metadata_df, index_column='File.Name', key_column='Corrupted', key_value=True):
    """
    Indexes metadata on specified column (index_column), then filters rows from the input data DataFrame 
    based on bool values in a specific metadata column (key_column).

    If the specified index_column is not found in metadata, logs a ERROR: INDEX KEY ERROR, 
    and returns unfiltered data.

    If the specified key_column is not found after indexing metadata, logs a ERROR: COLUMN KEY ERROR, 
    and returns unfiltered data.

    author: SB, BZ

    Parameters
    ----------
    data_df : pd.DataFrame
        DataFrame containing the data, indexed by features/precursors, with columns as samples.
    metadata_df : pd.DataFrame
        DataFrame containing metadata, indexed by sample.
    index_column : str
        Column name in meta_df to set as index to align with data_df.
    key_column : str
        Column name in meta_df to check for corrupted samples due to MS acquisition operator 
        technicalities. Default is 'Corrupted'. 
    key_value : bool or str
        the value in key_column specifying samples to be dropped.

    Returns
    -------
    pd.DataFrame
        Filtered version of data_df where only samples that lack corruption in metadata are retained 
        (i.e. Corrupt != True).
    log : str
        log string, logged with ###ERROR### messages if index and column keys not found in metadata
    """
    
    if index_column not in metadata_df.columns:
        log = f"###ERROR###: INDEX KEY ERROR - index_column '{index_column}' not found in metadata. Returning unfiltered data.\n"
        return data_df, log

    metadata_df_indexed = metadata_df.set_index(index_column)

    if key_column not in metadata_df_indexed.columns:
        log = f"###WARNING###: COLUMN KEY WARNING - key_column '{key_column}' not found in metadata. Returning unfiltered data.\n"
        return data_df, log
    
    log = "" #f"{timestamp()}\tAll keys found ('{index_column}', '{key_column}') in metadata. Filtering wide by metadata.\n"
    # Define columns to be dropped by samples with key_value in key_column intersecting with data_df
    drop_cols = list(set(metadata_df_indexed[metadata_df_indexed[key_column]==key_value].index) & set(data_df.columns))
    data_df = data_df.drop(drop_cols,axis=1)
    #remove precursors with all missing values after sample exclusion
    data_df = data_df.dropna(how="all",axis=0)
    return data_df, log


def filter_samples_by_detection_proportion(data_df, sample_thresh=None, sample_alpha=0.95):
    """
    Integrated pipeline that filters samples from the dataset based on a calculated or 
    user-defined optimal missingness threshold.

    author: SB, BZ

    Parameters
    ----------
    data_df : pd.DataFrame
        Input DataFrame (indexed on precursors)
    sample_thresh : float or None, default None
        Missingness threshold to use, between 0 and 1; if no hard threshold specified, triggers
        default (sample_thresh = None) for automated thresholding using `calculate_optimum_threshold`.
    sample_alpha : float 
        Alpha (significance) level for threshold calculation between 0 and 1 (default 0.95).
    
    Returns
    -------
    data_df : pd.DataFrame
        Filtered DataFrame retaining samples with missingness below the selected threshold.
    log : str
       Log string with error and execution comments.
        
    """
    log = ""
    
    # Validate allowed values for sample_thresh
    if sample_thresh is not None:
        if not isinstance(sample_thresh, (int, float)):
            log += f"###ERROR###: Invalid type for sample_thresh. Expected int or float, got {type(sample_thresh).__name__}.\n"
            return data_df, log
        if not (0 <= sample_thresh <= 1):
            log += f"###ERROR###: Invalid value for sample_thresh. Must be between 0 and 1. Got {sample_thresh}.\n"
            return data_df, log
        
        missingness = condition_on_detection_proportion(data_df)
        
        optimal_threshold = sample_thresh * 100
        log += f"\t\t\t\tUsing user-defined thresholding at {optimal_threshold:.2f}% missingness for sample filtering.\n"

    else:
        # Validate allowed values for alpha_sample
        if not isinstance(sample_alpha, (int, float)):
            log += f"###ERROR###: Invalid type for sample_alpha. Expected int or float, got {type(sample_alpha).__name__}.\n"
            return data_df, log
        if not (0 <= sample_alpha <= 1):
            log += f"###ERROR###: Invalid value for sample_alpha. Must be between 0 and 1. Got {sample_alpha}.\n"
            return data_df, log
        
        #calculate optimal sample threshold
        missingness, optimal_threshold = calculate_optimum_threshold(data_df, sample_alpha)

        log += f"\t\t\t\tUsing automated thresholding {optimal_threshold:.2f}% missingness (alpha={sample_alpha}) for sample filtering.\n"
    
    #only keep samples with missingness below optimal_threshold
    data_df = data_df.loc[:, missingness < optimal_threshold]
    #remove precursors with all missing values after sample exclusion
    data_df = data_df.dropna(how="all",axis=0)
    return data_df, log

##############################################################################
# Precursor and Protein Thresholding Functions
##############################################################################

def filter_precursors_by_detection_proportion(data_df, id_mapping_df=None, precursor_thresh=None, precursor_alpha=0.95, filter_per_protein=False, protein_col="Protein.Group", precursor_col="Precursor.Id"):
    """
    Integrated pipeline that filters precursors from the dataset based on a calculated or 
    user-defined optimal missingness threshold.

    author: SB, BZ

    Parameters
    ----------
    data_df : pd.DataFrame
        Input DataFrame (indexed on precursors).
    id_mapping_df: pd.DataFrame
        dataframe that stores mapping of precursors and proteins. 
    precursor_thresh : float or None, default None
        Missingness threshold to use (0–100%). If None, uses automatic thresholding via 
        `calculate_optimum_threshold`.
    precursor_alpha : float
        Alpha (significance) level for threshold calculation between 0 and 1(default is 0.95).
    filter_per_protein : bool
        if True, apply threshold protein-wise only to samples containing any precursor of a protein, if False apply filtering considering always all samples
    protein_col : str
        Name of the ID column containing the protein or gene ID. The default is "Protein.Group".
    precursor_col : str
        Name of the ID column in id_mapping_df, containing the precursor ID. The default is "Precursor.Id".
        
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame retaining precursors with missingness below the selected threshold.
    log : str
        Log string to pick up run time comments and warnings.
    """
    
    #logging
    log=""
    
    if filter_per_protein and id_mapping_df is None:
        log += "###ERROR###: No id_mapping_df passed while filter_per_protein set to True. Either set filter_per_protein to False for global thresholding or provide an id_mapping_df to map precursors to protein groups.\n"
    
    # Validate allowed values for precursor_thresh
    if precursor_thresh is not None:
        if not isinstance(precursor_thresh, (int, float)):
            log += f"###ERROR###: Invalid type for precursor_thresh. Expected int or float, got {type(precursor_thresh).__name__}.\n"
            return data_df, log
        if not (0 <= precursor_thresh <= 1):
            log += f"###ERROR###: Invalid value for precursor_thresh. Must be between 0 and 1. Got {precursor_thresh}.\n"
            return data_df, log
        
        #Calculate missingness
        missingness = condition_on_detection_proportion(data_df.T)
        optimal_threshold = precursor_thresh*100
        log += f"\t\t\t\tUsing user-defined thresholding at {optimal_threshold:.2f}% missingness for precursor filtering.\n"
        
    else:
        # Validate allowed values for alpha_pr
        if not isinstance(precursor_alpha, (int, float)):
            log += f"###ERROR###: Invalid type for precursor_alpha. Expected int or float, got {type(precursor_alpha).__name__}.\n"
            return data_df, log
        if not (0 <= precursor_alpha <= 1):
            log += f"###ERROR###: Invalid value for precursor_alpha. Must be between 0 and 1. Got {precursor_alpha}.\n"
            return data_df, log
    
        # Calculate missingness and global optimal threshold
        missingness, optimal_threshold = calculate_optimum_threshold(data_df.T, precursor_alpha, statistics='non-parametric')
        log += f"\t\t\t\tUsing automated thresholding {optimal_threshold:.2f}% missingness (alpha={precursor_alpha}) for precursor filtering.\n"
        
    # Filter by missingness
    if filter_per_protein:
        log += f"\t\t\t\tFiltering by missingness per protein according to {protein_col} across samples containing at least one precursor of the same protein.\n"
        keep_precursors = list()
        for pg in id_mapping_df[id_mapping_df[precursor_col].isin(data_df.index)][protein_col].unique():
            precursors_pg = list(set(id_mapping_df[id_mapping_df[protein_col]==pg][precursor_col]) & set(data_df.index))
            pg_threshold = 100 - (100-optimal_threshold) * (sum(data_df.loc[precursors_pg].sum()>0) / data_df.shape[1])
            keep_precs_df = missingness.loc[precursors_pg] < pg_threshold
            keep_precursors += keep_precs_df[keep_precs_df].index.to_list()
        data_df = data_df.loc[keep_precursors]
    else:
        data_df = data_df.loc[missingness < optimal_threshold]
        
    return data_df, log


def filter_precursors_by_exclusionlist(data_df,id_mapping_df,id_col,exclusion_id_list,precursor_col="Precursor.Id"):
    """
    Takes a list with IDs and excludes the corresponding precursors from the data_df.
    
    author: BZ
    
    Parameters
    ----------
    data_df : pd.DataFrame
        Input DataFrame (indexed on precursors).
    id_mapping_df: pd.DataFrame
        dataframe that stores mapping of precursors and proteins. 
    id_col : str
        column of id_mapping_df containing the ids from exclusion_id_list.
    exclusion_id_list : iterable
        iterable containing the ids to be excluded from data_df.
    precursor_col : str, optional
        Name of the column in id_mapping_df matching data_df.index. The default is "Precursor.Id".

    Returns
    -------
    pd.DataFrame
        filtered dataframe without ids provided by exclusion_id_list.
    log: str
        Log string to pick up run time comments and warnings.

    """
    log = ""
    if id_col not in id_mapping_df.columns:
        log += f"###ERROR###: COLUMN KEY ERROR - column '{id_col}' not found in id_mapping_df. Returning unfiltered data.\n"
    if precursor_col not in id_mapping_df.columns:
        log += f"###ERROR###: COLUMN KEY ERROR - column '{precursor_col}' not found in id_mapping_df. Returning unfiltered data.\n"
        
    # catch error and return
    if "###ERRROR###" in log:
        return data_df, log
    else:
        exclude_precursors = list(set(data_df.index) & set(id_mapping_df[id_mapping_df[id_col].isin(exclusion_id_list)][precursor_col]))
        exclusion_ids_in_data = id_mapping_df.loc[id_mapping_df[precursor_col].isin(exclude_precursors),id_col].unique()
        log += f"\t\t\t\t{len(exclusion_ids_in_data)} {id_col} out of {len(exclusion_id_list)} in exclusion list filtered from data matrix.\n"
        return data_df.drop(exclude_precursors,axis=0), log


def filter_proteins_by_id_count(data_df, id_mapping_df, min_ids=1, protein_col="Protein.Group", precursor_col="Precursor.Id", id_col="Peptide.Id"):
    """
    Remove proteins identified by less ids (precursors or peptides) than specified with min_ids.
    E.g. For filtering on a minimum number of identified precursors, set id_col="Precursor.Id" (same as precursor_col).
    E.g. For filtering on a minimum number of identified peptides, set id_col="Stripped.Sequence" (a different column in id_mapping_df).
    Note: Stripped.Sequence only serves as a proxy for peptides if no variable modifications are included in the dataset, i.e. a precursor can only correspond to a peptide with a specific charge state. 

    author: BZ, SB
    
    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    id_mapping_df: pd.DataFrame
        dataframe that stores mapping of precursors and proteins. 
    min_ids : int
        minimum number of ids required per protein, proteins with less ids will be excluded. The default is 1. 
    protein_col : str
        Name of the ID column containing the protein or gene ID. The default is "Protein.Group".
    precursor_col : str
        Name of the ID column in id_mapping_df, containing the precursor ID. The default is "Precursor.Id".
    id_col : str
        Name of the ID column in id_mapping_df, containing the ID used for counting. The default is "Peptide.Id".

    Returns
    -------
    data_df : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    log : str
        Log string to pick up run time comments and warnings.
    """
    
    log = ""
    
    # Check whether defined columns are present
    if id_col not in id_mapping_df.columns:
        log += f"##ERROR###: COLUMN KEY ERROR - column '{id_col}' not found in id_mapping_df. Returning unfiltered data.\n"
    if precursor_col not in id_mapping_df.columns:
        log += f"###ERROR###: COLUMN KEY ERROR - column '{precursor_col}' not found in id_mapping_df. Returning unfiltered data.\n"
    if protein_col not in id_mapping_df.columns:
        log += f"###ERROR###: COLUMN KEY ERROR - column '{precursor_col}' not found in id_mapping_df. Returning unfiltered data.\n"
    
    # catch error and return
    if "###ERRROR###" in log:
        return data_df, log
        
    if id_col==precursor_col: # filter by number of precursors per protein
        # Count precursors / protein present in data_df
        protein_counts = id_mapping_df[id_mapping_df[precursor_col].isin(data_df.index)][protein_col].value_counts()
    else: # filter by number of peptides per protein
        # Count precursors / id present in data_df
        precursor_counts = id_mapping_df[id_mapping_df[precursor_col].isin(data_df.index)].groupby(by=[id_col])[protein_col].value_counts()
        # Count ids / protein
        protein_counts = precursor_counts.index.get_level_values(protein_col).value_counts()
    
    # Keep proteins with at least n >= min_peptides
    proteins_to_keep = protein_counts[protein_counts >= min_ids].index

    # Get precursors mapped to filtered proteins
    precursors_to_keep = id_mapping_df[id_mapping_df[protein_col].isin(proteins_to_keep)][precursor_col].drop_duplicates()
    precursors_to_keep = precursors_to_keep[precursors_to_keep.isin(data_df.index)]

    # Filter on proteome level
    filtered_df = data_df.loc[data_df.index.isin(precursors_to_keep)]

    # Update logfile with dropped precursor count
    original_precursor_count = data_df.shape[0]
    filtered_precursor_count = filtered_df.shape[0]
    precursors_dropped = original_precursor_count - filtered_precursor_count
    original_protein_count = len(id_mapping_df[id_mapping_df[precursor_col].isin(data_df.index)][protein_col].unique())
    #filtered_protein_count = len(set(proteins_to_keep) & set(id_mapping_df[id_mapping_df[precursor_col].isin(data_df.index)][protein_col]))
    filtered_protein_count = len(id_mapping_df[id_mapping_df[precursor_col].isin(filtered_df.index)][protein_col].unique())
    proteins_dropped = original_protein_count - filtered_protein_count
    log += f"\t\t\t\tFiltered out {int(precursors_dropped)} precursors (kept {int(filtered_precursor_count)} out of {int(original_precursor_count)}) based on minimum of {min_ids} {id_col}s per protein.\n"
    log += f"\t\t\t\tFiltered out {int(proteins_dropped)} proteins (kept {int(filtered_protein_count)} out of {int(original_protein_count)}) based on minimum of {min_ids} {id_col}s per protein.\n"

    return filtered_df, log


#DEPREATED (new function filter_proteins_by_id_count to apply filtering on number of precursors or number of peptides)
'''
def filter_proteins_by_precursor_count(data_df, id_mapping_df, min_precursors=1, protein_col="Protein.Group", precursor_col="Precursor.Id"):
    """
    Remove proteins identified by less precursors than specified with min_precursors.

    author: BZ, SB
    
    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    id_mapping_df: pd.DataFrame
        dataframe that stores mapping of precursors and proteins. 
    min_precursors : int
        minimum number of precursors required per protein, proteins with less precursors will be excluded. The default is 1. 
    protein_col : str
        Name of the ID column containing the protein or gene ID. The default is "Protein.Group".
    precursor_col : str
        Name of the ID column in id_mapping_df, containing the precursor ID. The default is "Precursor.Id".

    Returns
    -------
    data_df : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    log : str
        Log string to pick up run time comments and warnings.
    """
    # Count precursors / protein present in data_df
    protein_counts = id_mapping_df[id_mapping_df[precursor_col].isin(data_df.index)][protein_col].value_counts()

    # Keep proteins with at least n >= min_precursors
    proteins_to_keep = protein_counts[protein_counts >= min_precursors].index

    # Get precursors mapped to filtered proteins
    precursors_to_keep = id_mapping_df[id_mapping_df[protein_col].isin(proteins_to_keep)][precursor_col].drop_duplicates()
    precursors_to_keep = precursors_to_keep[precursors_to_keep.isin(data_df.index)]

    # Filter on proteome level
    filtered_df = data_df.loc[data_df.index.isin(precursors_to_keep)]

    # Update logfile with dropped precursor count
    timestamp = lambda: time.strftime("%Y-%m-%d%t%H:%M:%S")
    original_precursor_count = data_df.shape[0]
    filtered_precursor_count = filtered_df.shape[0]
    precursors_dropped = original_precursor_count - filtered_precursor_count
    original_protein_count = len(id_mapping_df[id_mapping_df[precursor_col].isin(data_df.index)][protein_col].unique())
    #filtered_protein_count = len(set(proteins_to_keep) & set(id_mapping_df[id_mapping_df[precursor_col].isin(data_df.index)][protein_col]))
    filtered_protein_count = len(id_mapping_df[id_mapping_df[precursor_col].isin(filtered_df.index)][protein_col].unique())
    proteins_dropped = original_protein_count - filtered_protein_count
    log = f"{timestamp()}\tFiltered out {int(precursors_dropped)} precursors (kept {int(filtered_precursor_count)} out of {int(original_precursor_count)}) based on minimum of {min_precursors} precursors per protein.\n"
    log += f"\t\t\t\tFiltered out {int(proteins_dropped)} proteins (kept {int(filtered_protein_count)} out of {int(original_protein_count)}) based on minimum of {min_precursors} precursors per protein.\n"

    return filtered_df, log
'''

#DEPRECATED
'''
def filter_precursors_by_metadata(data_df, id_mapping_df, precursor_col="Precursor.Id", filter_col="Proteotypic", filter_value=True):
    """
    Exclude precursors in data_df with filter_col not equalling filter_value in id_mapping_df (precursor metadata). Can be used to filter i.e. all non-proteotypic precursors

    Parameters
    ----------
    data_df : pd.DataFrame
        wide-format profiles, columns are samples
    id_mapping_df : pd.DataFrame
        dataframe that stores mapping of precursors and proteins and additional metadata for precursors such as i.e. proteotypicity. 
    precursor_col : str, optional
        Name of the ID column in id_mapping_df, containing the precursor ID. The default is "Precursor.Id".
    filter_col : str, optional
        name of column in id_mapping_df to perform filter on. The default is "Proteotypic".
    filter_value : str/int/float/..., optional
        value in filter_col specifying precursors to keep.

    Returns
    -------
    filtered_df : pd.DataFrame
        filtered wide-format profiles, columns are samples.

    """
    log = ""
    
    #check if filter column is present
    if filter_col not in id_mapping_df.columns:
        log += f"###ERROR###: COLUMN KEY ERROR - column '{filter_col}' not found in id_mapping_df. Returning unfiltered data.\n"
    if precursor_col not in id_mapping_df.columns:
        log += f"###ERROR###: COLUMN KEY ERROR - column '{precursor_col}' not found in id_mapping_df. Returning unfiltered data.\n"
    
    # catch error and return
    if "###ERRROR###" in log:
        return data_df, log
    
    #apply filtering
    keep_precursors = list(set(data_df.index) & set(id_mapping_df[id_mapping_df[filter_col]==filter_value][precursor_col]))
    return([data_df.loc[keep_precursors], log])
'''

def filter_precursors_by_metadata(data_df, id_mapping_df, precursor_col="Precursor.Id", filter_col="Proteotypic", filter_value=True, operator="=="):
    """
    Exclude precursors in data_df with filter_col not equalling filter_value in id_mapping_df (precursor metadata). Can be used to filter i.e. all non-proteotypic precursors

    Parameters
    ----------
    data_df : pd.DataFrame
        wide-format profiles, columns are samples
    id_mapping_df : pd.DataFrame
        dataframe that stores mapping of precursors and proteins and additional metadata for precursors such as i.e. proteotypicity. 
    precursor_col : str, optional
        Name of the ID column in id_mapping_df, containing the precursor ID. The default is "Precursor.Id".
    filter_col : str, optional
        name of column in id_mapping_df to perform filter on. The default is "Proteotypic".
    filter_value : str/int/float/..., optional
        value in filter_col specifying precursors to keep.
    operator : str, optional
        string that defines the operator for the filtering, e.g. "==" to filter for filter_col == filter_value or "!=" to filter for filter_col != filter_value

    Returns
    -------
    filtered_df : pd.DataFrame
        filtered wide-format profiles, columns are samples.

    """
    log = ""
    
    #check correct data types before running eval() for filtering
    if type(id_mapping_df) != pd.DataFrame:
        log += "###ERROR###: id_mapping_df provided to filter_precursors_by_metadata is not a dataframe.\n"
    if operator not in ["==","!=","<",">","<=",">=","&","|"]:
        log += "###ERROR###: Invalid operator provided to filter_precursors_by_metadata, please use one of ==, !=, <, >, <=, >=, &, | .\n"
    #check if filter column is present
    if filter_col not in id_mapping_df.columns:
        log += f"###ERROR###: COLUMN KEY ERROR - column '{filter_col}' not found in id_mapping_df. Returning unfiltered data.\n"
    if precursor_col not in id_mapping_df.columns:
        log += f"###ERROR###: COLUMN KEY ERROR - column '{precursor_col}' not found in id_mapping_df. Returning unfiltered data.\n"
    
    # catch error and return
    if "###ERRROR###" in log:
        return data_df, log
    
    #apply filtering
    try:
        keep_precursors = list(set(data_df.index) & set(eval("id_mapping_df[id_mapping_df[filter_col]"+operator+str(filter_value)+"][precursor_col]")))
    except:
        try:
            keep_precursors = list(set(data_df.index) & set(eval("id_mapping_df[id_mapping_df[filter_col]"+operator+repr(filter_value)+"][precursor_col]")))
        except:
            log += "###ERROR###: Filtering on metadata failed for "+f"id_mapping_df[id_mapping_df[{filter_col}]{operator}{str(filter_value)}]{[precursor_col]}"+".\n"
            log += "###ERROR###: Filtering on metadata failed for "+f"id_mapping_df[id_mapping_df[{filter_col}]{operator}{repr(filter_value)}]{[precursor_col]}"+".\n"
    return([data_df.loc[keep_precursors], log])


##############################################################################
#  Precursor label filtering and replacement
##############################################################################

def filter_precursors_by_peptide_labels(data_df, remove_label=None, regex=True):
    """
    Specify a e.g. SILAC label (part of the Precursor name), all precursors containing the label will be removed from further analysis

    author: BZ, SB
    
    Parameters
    ----------
    data_df : pd.DataFrame
        wide-format profiles, columns are samples
    remove_label : str
        name of the label to identify precursors to be removed, default = None
    regex : bool
        if True, remove_label will be interpreted as a regular expression, if False as literal string. The default is True.

    Returns
    -------
    pd.DataFrame
        wide-format profiles, columns are samples, index with Precursor IDs.
    log : str
        Log string to pick up run time comments and warnings.
    """
    
    log = ""
    
    if not remove_label:
        log += "###WARNING###: No precursors removed: remove_label is empty or None.\n"
        return data_df, log
    
    #obtain list of precursors to be excluded
    flagged_precursors = data_df.index[data_df.index.str.contains(remove_label,regex=regex)]
    
    if regex:
        appendix = " interpreted as regular expression"
    else:
        appendix = ""
        
    if len(flagged_precursors)>0:
        log += f"\t\t\t\t{len(flagged_precursors)} precursors excluded because they contained the label: '{remove_label}'{appendix}.\n"
    else:
        log += f"\t\t\t\tNo precursors excluded: label '{remove_label}'{appendix} not found in any index.\n"
    
    return data_df.drop(flagged_precursors, axis=0), log



def replace_labelled_precursors(data_df, id_mapping_df, mask_df, precursor_col="Precursor.Id", target_label=None, replace_label=None): 
    """
    Specify a e.g. SILAC label (part of the Precursor name), all precursors containing the label will be replaced with another label.

    author: BZ, SB
    
    Parameters
    ----------
    data_df : pd.DataFrame
        wide-format profiles, columns are samples
    id_mapping_df: pd.DataFrame
        dataframe that stores mapping of precursors and proteins. 
    mask_df : pd.DataFrame
        DataFrame storing mask values for features or samples.
    precursor_col : str
        Name of the ID column in id_mapping_df, containing the precursor ID. The default is "Precursor.Id".
    target_label : str
        name of the label to identify precursors to be replaced. Default is None.
    replace_label : str
        Substring to replace the target_label with. Default is None.

    Returns
    -------
    pd.DataFrame
        wide-format profiles, columns are samples, index with Precursor IDs.
    log : str
        Log string to pick up run time comments and warnings.
    """
    
    log = ""
    
    if not target_label:
        log += "###WARNING###: No replacement performed: target_label is empty or None.\n"
    if replace_label == None:
        log += "###WARNING###: No replacement performed: replace_label is None.\n"  
    if "###WARNING###" in log:
        return data_df, log
    
    # collect all precursors across data, mask and id_mapping_df
    precursors = set(data_df.index) | set(mask_df.index) | set(id_mapping_df[precursor_col])
    
    # Find matching precursors
    flagged_precursors = list(filter(lambda x: target_label in x, precursors))
    
    if len(flagged_precursors)>0:
        # Create a dataframe for translations
        replace_df = pd.DataFrame(index=flagged_precursors,columns=["replaced_precursor"])
        replace_df["replaced_precursor"] = replace_df.index.to_series().replace(to_replace=target_label, value=replace_label, regex=True)
        # Update indices
        data_df.rename(index=replace_df["replaced_precursor"].to_dict(),inplace=True)
        mask_df.rename(index=replace_df["replaced_precursor"].to_dict(),inplace=True)
        id_mapping_df[precursor_col] = id_mapping_df[precursor_col].replace(replace_df["replaced_precursor"].to_dict())
        # Update log
        num_replaced_precs_in_data = len(set(data_df.index) & set(replace_df["replaced_precursor"]))
        log += f"\t\t\t\t{len(flagged_precursors)} precursors replaced: '{target_label}' -> '{replace_label}'\n"
        log += f"\t\t\t\tOf these {num_replaced_precs_in_data} precursors were present in the data matrix.\n"
    else:
        log += f"\t\t\t\tNo replacement performed: {target_label} was not found in any of the precursors in data matrix, mask nor id mapping.\n"

    return data_df, log




##############################################################################
#  Other Functions (Automated Thresholding, Mask Updater)
##############################################################################

def condition_on_detection_proportion(data_df):

    """
    Calculates the percentage of missing values for each sample in the input DataFrame.

    author: SB

    Parameters
    ----------
    data_df : pd.DataFrame
        Wide format DataFrame indexed by features/precursors
    
    Returns
    -------
    missingness : pd.Series
        Series indexed by samples, indicating percentage of missing values.
    """
    missing = data_df.isna().mean(axis=0) * 100
    missingness = missing.sort_values()
    return missingness


def calculate_optimum_threshold(data_df, alpha=0.95, statistics='parametric'):
    """
    Calculates the optimal missingness threshold based on a statistical approach, 
    parametric (normal distribution), by default.

    author: SB

    Parameters
    ----------
    input_df : pd.DataFrame
        Input DataFrame containing missing values.
    alpha : float, optional
        Confidence interval alpha level for threshold calculation (default 0.95).

    Returns
    -------
    missingness : pd.Series
        Series indicating percentage of missing values per precursor.
    upper_bound : float
        Calculated threshold value (upper bound) for missingness percentage.
    """

    missingness = condition_on_detection_proportion(data_df)

    if statistics == 'parametric':
        mu = np.mean(missingness)
        sigma = np.std(missingness, ddof=1)
        z = stats.norm.ppf(alpha)  
        upper_bound = mu + z * sigma

    elif statistics == 'non-parametric':
        upper_bound = missingness.quantile(alpha)

    return missingness, upper_bound




def update_mask(data_df, axis, method, mask_df, mask_description_df, pre_method_columns=[], pre_method_rows=[]):
    """
    Update masking information and related metadata for filtered data.

    This function updates `mask_df`, `mask_description_df`, `metadata_df`, 
    and `id_mapping_df` based on differences between the original (`data_df`) 
    and processed (`data_processed_df`) datasets. It assigns a new mask value 
    for the given filtering `method` and updates relevant metadata depending 
    on the specified `axis`.

    The new mask value is determined as one greater than the 
    second-largest value in `mask_description_df['value']` 
    (or the largest if fewer than two values exist).
    
    author: SB, BZ

    Parameters
    ----------
    data_df : pd.DataFrame
        The processed dataset after filtering or transformation.
    axis : str
        Axis along which to apply masking. Must be one of:
        - `'samples'` : Compare columns between `data_df` and `data_processed_df`.
        - `'precursors'` : Compare index between `data_df` and `data_processed_df`.
    method : str
        The name of the filtering method being applied.
    mask_df : pd.DataFrame
        DataFrame storing mask values for features or samples.
    mask_description_df : pd.DataFrame
        DataFrame mapping filter methods to mask values. 
        Must contain columns: `['filter', 'value']`.

    Returns
    -------
    tuple of pd.DataFrame
        Updated `(mask_df, mask_description_df)`.

    Notes
    -----
    - When `axis='samples'`, masking is applied to columns present in `data_df` 
      but missing from `data_processed_df`.
    - When `axis='precursors'`, masking is applied to rows present in `data_df` 
      but missing from `data_processed_df`.
    """
    
    if method not in mask_description_df.index:
        # Get the second largest value in mask_description_df safely
        sorted_values = mask_description_df['value'].sort_values()
        if len(sorted_values) >= 2:
            second_max_mapping = sorted_values.iloc[-2]
        else:
            second_max_mapping = sorted_values.iloc[-1]
        new_value = second_max_mapping + 1
        mask_description_df.loc[method] = {'value': new_value}
    value = mask_description_df.loc[method,"value"]

    # Update based on axis
    if axis == 'samples':
        flags = pre_method_columns.difference(data_df.columns)
        if not flags.empty:
            mask_df[flags] = mask_df[flags].replace(mask_description_df.loc["measured and kept","value"],value)

    elif axis == 'precursors':
        flags = pre_method_rows.difference(data_df.index)
        if not flags.empty:
            mask_df.loc[flags,:] = mask_df.loc[flags,:].replace(mask_description_df.loc["measured and kept","value"],value)

    return mask_df, mask_description_df, flags



##############################################################################
#  Main
##############################################################################

if __name__ == "__main__":

    module = "FILTERING_WIDE"
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
    
    #load 
    data_df, parameters_dict, id_mapping_df, seed, log, metadata_df, mask_df, mask_description_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = module, 
                                                                                                                                                               get_metadata=True, 
                                                                                                                                                               get_id_mapping=True, 
                                                                                                                                                               get_mask=True, 
                                                                                                                                                               return_data_fp_only=False,
                                                                                                                                                               skip_parser=skip_parser,
                                                                                                                                                               data_df_fp=data_df_fp,
                                                                                                                                                               param_fp=param_fp,
                                                                                                                                                               out_dir_fp=out_dir_fp,
                                                                                                                                                               meta_df_fp=meta_df_fp,
                                                                                                                                                               log_file=log_file,
                                                                                                                                                               idmap_fp=idmap_fp,
                                                                                                                                                               counter=counter,
                                                                                                                                                               outname=outname)
 

    try:
        # Check if loading was successful, else pass error to log and exit   
        helper.logger(log_file,log)
        if metadata_df is None:
            helper.logger(log_file,"###ERROR###: Metadata loading failed.\n")
            if exit_module:
                sys.exit("###ERROR###: Metadata loading failed.")
        if "###WARNING###" in log:
            helper.logger(log_file,"###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"]+"\n")
            if exit_module:
                sys.exit("###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])
        
        #helper functions
        timestamp = lambda: time.strftime("%Y-%m-%d%t%H:%M:%S")
        
        # Log parameters
        #log = f"{timestamp()}\tFiltering Wide Format\n"
        #log += "\nParameters\n------------\n"
        #log += "\t" + pd.DataFrame(parameters_dict.items(),columns=["Parameter","Value"]).to_string(index=False).replace("\n","\n\t") 
        #log += "\n------------\n"  
        #helper.logger(log_file, module_output=log)
        
        #parameter extraction
        methods = parameters_dict['filtering_wide_methods'].replace(" ", "").split(";")
        #for multiple calls of filter_precursors_by_metadata define a counter accessing correct parameter subset, initialized with zero
        filter_precursors_by_metadata_counter = 0
        data_processed_df = data_df
        
        #create filtering metadata to remember excluded samples, excluded precursors will be flagged in id_mapping_df
        sample_filters_df = pd.DataFrame(index=data_processed_df.columns)
        sample_filters_df.index.name = parameters_dict['metadata_sample_column']
        
        if methods == []:
            log = "###WARNING###: No filtering methods provided. Skipping filtering on wide format.\n"
            helper.logger(log_file,log)
        else:
            for method in methods:
                
                # snapshot BEFORE this method
                pre_method_columns = data_processed_df.columns
                pre_method_rows = data_processed_df.index
                
                # SAMPLE FILTERS
                if method == "filter_samples_based_on_metadata":
                    #for every column specified for exclusion by metadata
                    for k,key_column in enumerate(parameters_dict['metadata_exclusion_columns'].split(";")):
                        # snapshot BEFORE this method
                        pre_method_columns = data_processed_df.columns
                        pre_method_rows = data_processed_df.index
                        
                        #select corresponding value for exclusion in the given key_column
                        key_value = eval(parameters_dict["metadata_exclusion_values_types"].split(";")[k])(str(parameters_dict['metadata_exclusion_values']).split(";")[k])
                        
                        #log filtering method
                        log = f"{timestamp()}\tExcluding samples from data with {key_value} in metadata column {key_column}.\n"
                        helper.logger(log_file,log)
                        #apply filtering based on metadata
                        data_processed_df, log  = filter_samples_based_on_metadata(data_df=data_processed_df, 
                                                                                   metadata_df=metadata_df, 
                                                                                   index_column=parameters_dict['metadata_sample_column'], 
                                                                                   key_column=key_column,
                                                                                   key_value=key_value)
                        #catch errors
                        if "###ERROR###" in log:
                            helper.logger(log_file,log)
                            if exit_module:
                                sys.exit(log)
                        if not "###WARNING###" in log:
                            # update mask_df
                            mask_df, mask_description_df, flags = update_mask(data_df=data_processed_df,
                                                                              axis='samples', 
                                                                              method=key_column+"="+str(key_value),
                                                                              mask_df=mask_df, 
                                                                              mask_description_df=mask_description_df, 
                                                                              pre_method_columns=pre_method_columns,
                                                                              pre_method_rows=pre_method_rows)
                            
                            #remember excluded samples
                            sample_filters_df.loc[flags,key_column+"="+str(key_value)] = True
                            
                            #log filtering progress
                            num_excluded_samples = len(pre_method_columns) - len(data_processed_df.columns)
                            log += f"\t\t\t\t{num_excluded_samples} samples excluded by {key_value} in metadata_column {key_column}.\n"
                        helper.logger(log_file,log)
                        
                elif method == "filter_samples_by_detection_proportion":
                    #log filtering method
                    log = f"{timestamp()}\tFiltering samples from data by missingness threshold.\n"
                    helper.logger(log_file,log)
                    #apply filtering of samples by detecction proportion (missingness)
                    data_processed_df, log = filter_samples_by_detection_proportion(data_df=data_processed_df, 
                                                                                    sample_thresh=parameters_dict['sample_missingness_thresh'], 
                                                                                    sample_alpha=parameters_dict['sample_missingness_alpha'])
                    #catch errors
                    if "###ERROR###" in log:
                        helper.logger(log_file,log)
                        if exit_module:
                            sys.exit(log)
                        
                    # update mask_df
                    mask_df, mask_description_df, flags = update_mask(data_df=data_processed_df,
                                                                      axis='samples', 
                                                                      method='sample missingness',
                                                                      mask_df=mask_df, 
                                                                      mask_description_df=mask_description_df, 
                                                                      pre_method_columns=pre_method_columns,
                                                                      pre_method_rows=pre_method_rows)
                    
                    #remember excluded samples
                    sample_filters_df.loc[flags,"sample missingness"] = True
                    
                    #log filtering progress
                    num_excluded_samples = len(pre_method_columns) - len(data_processed_df.columns)
                    percent_excluded_samples = np.round(num_excluded_samples/len(data_processed_df.columns)*100,2)
                    log += f"\t\t\t\t{num_excluded_samples} ({percent_excluded_samples}%) samples filtered due to too many missing values.\n"
                    helper.logger(log_file,log)
                
                elif method == "filter_by_normalization_factors":
                    pass #TODO
                    
                    
                # PRECURSOR FILTERS
                elif method == "filter_precursors_by_detection_proportion":
                    #log filtering method
                    log = f"{timestamp()}\tFiltering precursors from data by missingness threshold.\n"
                    helper.logger(log_file,log)
    
                    #apply filtering of precursors by detection proportion (missingness)
                    data_processed_df, log = filter_precursors_by_detection_proportion(data_df=data_processed_df, 
                                                                                       id_mapping_df=id_mapping_df,
                                                                                       precursor_thresh=parameters_dict["precursor_missingness_thresh"], 
                                                                                       precursor_alpha=parameters_dict["precursor_missingness_alpha"],
                                                                                       filter_per_protein=parameters_dict["proteinwise_precursor_missingness"],
                                                                                       protein_col=parameters_dict["protein_col"],
                                                                                       precursor_col=parameters_dict["precursor_col"])
                    #catch errors
                    if "###ERROR###" in log:
                        helper.logger(log_file,log)
                        if exit_module:
                            sys.exit(log)
                    
                    # update mask_df
                    mask_df, mask_description_df, flags = update_mask(data_df=data_processed_df,
                                                                      axis='precursors', 
                                                                      method='precursor missingness',
                                                                      mask_df=mask_df, 
                                                                      mask_description_df=mask_description_df, 
                                                                      pre_method_columns=pre_method_columns,
                                                                      pre_method_rows=pre_method_rows)
                    
                    #remember excluded precursors in id_mapping_df
                    id_mapping_df.loc[id_mapping_df['Precursor.Id'].isin(flags), 'precursor missingness'] = True
                    id_mapping_df['precursor missingness'] = id_mapping_df['precursor missingness'].fillna(False)
    
                    #log filtering progress
                    num_excluded_precursors = len(pre_method_rows) - len(data_processed_df.index)
                    percent_excluded_precursors = np.round(num_excluded_precursors/len(pre_method_rows)*100,2)
                    log += f"\t\t\t\t{num_excluded_precursors} ({percent_excluded_precursors}%) precursors excluded by exceeding missingness threshold.\n"
                    helper.logger(log_file,log)
    
                elif method == "filter_proteins_by_peptide_count":
                    #log filtering method
                    log = f"{timestamp()}\tFiltering proteins from data with less than {parameters_dict['min_peptides']} peptides per protein.\n"
                    helper.logger(log_file,log)
                    #apply filtering by minimum number of precursors per protein
                    data_processed_df, log = filter_proteins_by_id_count(data_df=data_processed_df, 
                                                                                id_mapping_df=id_mapping_df, 
                                                                                min_ids=parameters_dict["min_peptides"], 
                                                                                protein_col=parameters_dict["protein_col"], 
                                                                                precursor_col=parameters_dict["precursor_col"],
                                                                                id_col=parameters_dict["peptide_col"])
                    
                    #catch errors
                    if "###ERROR###" in log:
                        helper.logger(log_file,log)
                        if exit_module:
                            sys.exit(log)
    
                    #update mask_description_df
                    mask_description_df.index = mask_description_df.index.to_series().replace(to_replace= r"Min. .* Peptides", value="Min. " + str(parameters_dict["min_peptides"]) + " Peptides", regex=True)
    
                    # update mask_df
                    mask_df, mask_description_df, flags = update_mask(data_df=data_processed_df,
                                                                      axis='precursors', 
                                                                      method="Min. " + str(parameters_dict["min_peptides"]) + " Peptides",
                                                                      mask_df=mask_df, 
                                                                      mask_description_df=mask_description_df, 
                                                                      pre_method_columns=pre_method_columns,
                                                                      pre_method_rows=pre_method_rows)
                    
                    #remember excluded precursors in id_mapping_df
                    id_mapping_df.loc[id_mapping_df['Precursor.Id'].isin(flags), "Less than " + str(parameters_dict["min_peptides"]) + " Peptides"] = True
                    id_mapping_df["Less than " + str(parameters_dict["min_peptides"]) + " Peptides"] = id_mapping_df["Less than " + str(parameters_dict["min_peptides"]) + " Peptides"].fillna(False)
                    
                    #log filtering progress
                    helper.logger(log_file,log)
                    
                elif method == "filter_proteins_by_precursor_count":
                    #log filtering method
                    log = f"{timestamp()}\tFiltering proteins from data with less than {parameters_dict['min_precursors']} precursors per protein.\n"
                    helper.logger(log_file,log)
                    #apply filtering by minimum number of precursors per protein
                    data_processed_df, log = filter_proteins_by_id_count(data_df=data_processed_df, 
                                                                         id_mapping_df=id_mapping_df, 
                                                                         min_ids=parameters_dict["min_precursors"], 
                                                                         protein_col=parameters_dict["protein_col"], 
                                                                         precursor_col=parameters_dict["precursor_col"],
                                                                         id_col=parameters_dict["precursor_col"])
                    
                    #catch errors
                    if "###ERROR###" in log:
                        helper.logger(log_file,log)
                        if exit_module:
                            sys.exit(log)
    
                    #update mask_description_df
                    mask_description_df.index = mask_description_df.index.to_series().replace(to_replace= r"Min. .* Precursors", value="Min. " + str(parameters_dict["min_precursors"]) + " Precursors", regex=True)
    
                    # update mask_df
                    mask_df, mask_description_df, flags = update_mask(data_df=data_processed_df,
                                                                      axis='precursors', 
                                                                      method="Min. " + str(parameters_dict["min_precursors"]) + " Precursors",
                                                                      mask_df=mask_df, 
                                                                      mask_description_df=mask_description_df, 
                                                                      pre_method_columns=pre_method_columns,
                                                                      pre_method_rows=pre_method_rows)
                    
                    #remember excluded precursors in id_mapping_df
                    id_mapping_df.loc[id_mapping_df['Precursor.Id'].isin(flags), "Less than " + str(parameters_dict["min_precursors"]) + " Precursors"] = True
                    id_mapping_df["Less than " + str(parameters_dict["min_precursors"]) + " Precursors"] = id_mapping_df["Less than " + str(parameters_dict["min_precursors"]) + " Precursors"].fillna(False)
                    
                    #log filtering progress
                    helper.logger(log_file,log)
    
                elif method == "filter_precursors_by_metadata":
                    #log filtering method
                    log = f"{timestamp()}\tExcluding precursors from id_mapping_df with {parameters_dict['precursor_filter_col']} not {parameters_dict['precursor_filter_value']}.\n"
                    helper.logger(log_file,log)
                    #extract correct parameter subset
                    filter_col=parameters_dict["precursor_filter_col"].split(";")[filter_precursors_by_metadata_counter]
                    filter_value=parameters_dict["precursor_filter_value"].split(";")[filter_precursors_by_metadata_counter]
                    operator=parameters_dict["filter_operator"].split(";")[filter_precursors_by_metadata_counter]
                    #apply exclusion of selected precursors
                    data_processed_df, log = filter_precursors_by_metadata(data_df=data_processed_df, 
                                                                           id_mapping_df=id_mapping_df, 
                                                                           precursor_col=parameters_dict["precursor_col"], 
                                                                           filter_col=filter_col, 
                                                                           filter_value=filter_value,
                                                                           operator=operator)
                    
                    #catch errors
                    if "###ERROR###" in log:
                        helper.logger(log_file,log)
                        if exit_module:
                            sys.exit(log)
                        
                    # update mask_df
                    mask_df, mask_description_df, flags = update_mask(data_df=data_processed_df,
                                                                      axis='precursors', 
                                                                      method=parameters_dict["precursor_filter_col"],
                                                                      mask_df=mask_df, 
                                                                      mask_description_df=mask_description_df, 
                                                                      pre_method_columns=pre_method_columns,
                                                                      pre_method_rows=pre_method_rows)
                    
                    #log filtering progress
                    helper.logger(log_file,log)
                    
                elif method == "filter_precursors_by_exclusion_list":
                    #log filtering method
                    log = f"{timestamp()}\tFiltering {parameters_dict['exclusion_list_id_col']} from data as provided by {parameters_dict['exclusion_list_fp']}.\n"
                    helper.logger(log_file,log)
                    # Load exclusion list
                    if os.path.exists(parameters_dict["exclusion_list_fp"]):
                        id_list = pd.read_csv(parameters_dict["exclusion_list_fp"],header=None)[0].to_list()
                    #apply filtering by exclusion list
                    data_processed_df, log = filter_precursors_by_exclusionlist(data_df=data_processed_df,
                                                                                id_mapping_df=id_mapping_df,
                                                                                id_col=parameters_dict["exclusion_list_id_col"],
                                                                                exclusion_id_list=id_list,
                                                                                precursor_col=parameters_dict["precursor_col"])
                    #catch errors
                    if "###ERROR###" in log:
                        helper.logger(log_file,log)
                        if exit_module:
                            sys.exit(log)
                    
                    # update mask_df
                    mask_df, mask_description_df, flags = update_mask(data_df=data_processed_df,
                                                                      axis='precursors', 
                                                                      method="Exclusion by "+parameters_dict["exclusion_list_fp"],
                                                                      mask_df=mask_df, 
                                                                      mask_description_df=mask_description_df, 
                                                                      pre_method_columns=pre_method_columns,
                                                                      pre_method_rows=pre_method_rows)
                    
                    #remember excluded precursors in id_mapping_df
                    id_mapping_df.loc[id_mapping_df['Precursor.Id'].isin(flags), "Exclusion by "+parameters_dict["exclusion_list_fp"]] = True
                    id_mapping_df["Exclusion by "+parameters_dict["exclusion_list_fp"]] = id_mapping_df["Exclusion by "+parameters_dict["exclusion_list_fp"]].fillna(False)
                    
                    #log filtering progress
                    helper.logger(log_file,log)
                    
                elif method == "filter_precursors_by_peptide_labels":
                    #log filtering method
                    log = f"{timestamp()}\tFiltering precursors from data by defined precursor labels.\n"
                    helper.logger(log_file,log)
                    #apply filtering by pre-defined precursor labels
                    data_processed_df, log = filter_precursors_by_peptide_labels(data_df=data_processed_df, 
                                                                                 remove_label=parameters_dict["remove_label"], 
                                                                                 regex=parameters_dict["is_remove_label_regex"])
                    
                    # update mask_df
                    mask_df, mask_description_df, flags = update_mask(data_df=data_processed_df,
                                                                      axis='precursors', 
                                                                      method='precursor label: '+parameters_dict['remove_label'],
                                                                      mask_df=mask_df, 
                                                                      mask_description_df=mask_description_df, 
                                                                      pre_method_columns=pre_method_columns,
                                                                      pre_method_rows=pre_method_rows)
                    
                    #remember excluded precursors in id_mapping_df
                    id_mapping_df.loc[id_mapping_df['Precursor.Id'].isin(flags), 'precursor label: '+parameters_dict['remove_label']] = True
                    id_mapping_df['precursor label: '+parameters_dict['remove_label']] = id_mapping_df['precursor label: '+parameters_dict['remove_label']].fillna(False)
                    
                    #log filtering progress
                    helper.logger(log_file,log)
                
                # PRECURSOR RENAMING
                elif method == "replace_labelled_precursors":
                    #log
                    log = f"{timestamp()}\tReplacing precursors names in data with defined precursor labels.\n"
                    helper.logger(log_file,log)
                    
                    target_labels = parameters_dict["target_labels"].split(";")
                    replace_labels = parameters_dict["replace_labels"].split(";")
                    for i in range(len(target_labels)):
                        data_processed_df, log = replace_labelled_precursors(data_df=data_processed_df, 
                                                                             id_mapping_df=id_mapping_df,
                                                                             mask_df=mask_df,
                                                                             precursor_col=parameters_dict["precursor_col"],
                                                                             target_label=target_labels[i], 
                                                                             replace_label=replace_labels[i])
                        #log progress
                        helper.logger(log_file,log)
                    
                else:
                    log = f"###WARNING###: Unknown filtering method: {method} — skipping.\n"
                    helper.logger(log_file,log)
    
    
       
        
       
        #sort mask_description_df 
        mask_description_df.sort_values(by="value",inplace=True)
        module_output = ""
        #save sample_filters_df and fillna with False
        if not sample_filters_df.empty:
            sample_filters_df.fillna(False,inplace=True)
            sample_filters_df.to_csv(os.path.join(out_dir_fp, "filtered_samples.tsv"), sep="\t")
            module_output += f"{timestamp()}\tAdditional metadata containing information about filtered samples saved to {os.path.join(out_dir_fp, 'filtered_samples.tsv')}\n"
        #save id_mapping_df
        id_mapping_df.to_csv(os.path.join(out_dir_fp, "id_mapping.tsv"), sep="\t") # is this the correct location #???
        module_output += f"{timestamp()}\tInformation about filtered precursors saved to {os.path.join(out_dir_fp, 'id_mapping.tsv')}\n"
        
        #close module
        helper.close_module(data_processed_df=data_processed_df, masked_df=mask_df, mask_description_df=mask_description_df, log_file=log_file, module_output=module_output, out_dir_fp=out_dir_fp, outname=outname, mask=True, exit_module=exit_module)
    except Exception as e:
        e_type, e_obj, e_tb = sys.exc_info()
        script_name = os.path.split(e_tb.tb_frame.f_code.co_filename)[1]
        helper.logger(log_file, "###ERROR### Module " +script_name + " ran into an error at line : " + str(e_tb.tb_lineno) +"\n\t\t\t\t" + str(e) + "\n")
        if exit_module:
            sys.exit("Filtering_wide exited with an error.")
            
