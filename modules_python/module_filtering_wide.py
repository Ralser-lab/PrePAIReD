"""
Input: wide-format data
Output: Filtered wide-format data

functions:
    - filter_samples_based_on_metadata: Removes samples of poor quality, based on flags in a designated column in annotated metadata dataframe.
    - filter_samples_by_detection_proportion: Removes samples that contain missing values more than a specified or automated threshold value (% missing proportion).
    - filter_precursors_by_detection_proportion: Removes precursors that are detected less than a specified or automated threshold value (% missing proportion).
    - filter_precursors_by_peptide_labels: Removes precursors that contain specific substrings (e.g., SILAC labels like "(SILAC-K-L)").
    - filter_proteins_by_peptide_count: Removes precursors mapping to proteins less than a minimum number of unique peptides. 
    - condition_on_detection_proportion: Calculates % missing proportion by conditioning on missingness on either sample or precursor dimension.
    - calculate_optimum_threshold: Automatically determines the most appropriate detection threshold value (for precursor or sample by % missing proportion) based on sampling statistics.
    - replace_labelled_precursors: Finds and replaces specific labeling substrings (e.g., "(SILAC-K-L)") in precursor IDs with another label (e.g., "K"). 
    - update_mask_df: Updates mask_df with a unique identifier after a specific wide filtering step.
"""

import numpy as np
from datetime import datetime
import scipy.stats as stats
import module_helper as helper

##############################################################################
# Sample Thresholding Functions
##############################################################################

def filter_samples_based_on_metadata(data_df, meta_df, module_output, index_column='File.Name', key_column='Corrupted'):
    """
    Indexes metadata on specified column (index_column), then filters rows from the input data DataFrame 
    based on bool values in a specific metadata column (key_column).

    If the specified index_column is not found in metadata, logs a ERROR: INDEX KEY ERROR, 
    and returns unfiltered data.

    If the specified key_column is not found after indexing metadata, logs a ERROR: COLUMN KEY ERROR, 
    and returns unfiltered data.

    Author: SB

    Parameters
    ----------
    data_df : pd.DataFrame
        DataFrame containing the data, indexed by features/precursors, with columns as samples.
    meta_df : pd.DataFrame
        DataFrame containing metadata, indexed by sample.
    module_output : str
        Log file to pick up run time comments and warnings
    index_column : str
        Column name in meta_df to set as index to align with data_df.
    key_column : str
        Column name in meta_df to check for corrupted samples due to MS acquisition operator 
        technicalities. Default is 'Corrupted'. The column should be boolean.

    Returns
    -------
    pd.DataFrame
        Filtered version of data_df where only samples that lack corruption in metadata are retained 
        (i.e. Corrupt != True).
    str
        Module output, logged with ###ERROR### messages if index and column keys not found in metadata
    """
    timestamp = lambda: datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    if index_column not in meta_df.columns:
        module_output+=f"[{timestamp()}] ###ERROR###: INDEX KEY ERROR - index_column '{index_column}' not found in metadata. Returning unfiltered data.\n"
        return data_df, module_output

    meta_df_indexed = meta_df.set_index(index_column)

    if key_column not in meta_df_indexed.columns:
        module_output+=f"[{timestamp()}] ###ERROR###: COLUMN KEY ERROR - key_column '{key_column}' not found in metadata. Returning unfiltered data.\n"
        return data_df, module_output
    
    module_output+=f"[{timestamp()}] All keys found ('{index_column}', '{key_column}') in metadata. Filtering wide by metadata.\n"
    # Fill NAs in key column with FALSE
    sample_flags = meta_df_indexed.loc[:, key_column].fillna(False) 
    # Sample interesection between meta_df and data_df
    sample_flags = sample_flags[sample_flags.index.isin(data_df.columns)] 
    return data_df.loc[:, ~sample_flags], module_output

def filter_samples_by_detection_proportion(data_df, module_output, sample_thresh=None, sample_alpha=0.95):
    """
    Integrated pipeline that filters samples from the dataset based on a calculated or 
    user-defined optimal missingness threshold.

    author: SB

    Parameters
    ----------
    data_df : pd.DataFrame
        Input DataFrame (indexed on precursors)
    module_output : str
        Log file to pick up run time comments and warnings
    sample_thresh : float or None, default None
        Missingness threshold to use, between 0 and 100%; if no hard threshold specified, triggers
        default (optimum_sample = None) for automated thresholding using `calculate_optimum_threshold`.
    sample_alpha : float 
        Alpha (significance) level for threshold calculation between 0 and 1 (default 0.95).
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame retaining samples with missingness below the selected threshold.
    """
    timestamp = lambda: datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Validate allowed values for optimum_sample
    if sample_thresh is not None:
        if not isinstance(sample_thresh, (int, float)):
            module_output += f"[{timestamp()}] ###ERROR###: Invalid type for optimum_sample. Expected int or float, got {type(sample_thresh).__name__}.\n"
            return data_df, module_output
        if not (0 <= sample_thresh <= 100):
            module_output += f"[{timestamp()}] ###ERROR###: Invalid value for optimum_sample. Must be between 0 and 100. Got {sample_thresh}.\n"
            return data_df, module_output

    # Validate allowed values for alpha_sample
    if not isinstance(sample_alpha, (int, float)):
        module_output += f"[{timestamp()}] ###ERROR###: Invalid type for sample_alpha. Expected int or float, got {type(sample_alpha).__name__}.\n"
        return data_df, module_output
    if not (0 <= sample_alpha <= 1):
        module_output += f"[{timestamp()}] ###ERROR###: Invalid value for sample_alpha. Must be between 0 and 1. Got {sample_alpha}.\n"
        return data_df, module_output

    missingness, optimal_threshold = calculate_optimum_threshold(data_df, sample_alpha)

    if sample_thresh is not None:
        optimal_threshold = sample_thresh / 100 
        module_output += f"[{timestamp()}] Using hard thresholding at {sample_thresh} % missingness for sample filtering.\n"
    else:
        module_output += f"[{timestamp()}] Using automated thresholding {optimal_threshold:.2f}% missingness (alpha={sample_alpha}) for sample filtering.\n"

    return data_df.loc[:, missingness < optimal_threshold], module_output

##############################################################################
# Precursor and Protein Thresholding Functions
##############################################################################

def filter_precursors_by_detection_proportion(data_df, module_output, precursor_thresh=None, precursor_alpha=0.95):
    """
    Integrated pipeline that filters precursors from the dataset based on a calculated or 
    user-defined optimal missingness threshold.

    author: SB

    Parameters
    ----------
    data_df : pd.DataFrame
        Input DataFrame (indexed on precursors).
    module_output : str
        Log file to pick up run time comments and warnings
    precursor_thresh : float or None, default None
        Missingness threshold to use (0–100%). If None, uses automatic thresholding via 
        `calculate_optimum_threshold`.
    precursor_alpha : float
        Alpha (significance) level for threshold calculation between 0 and 1(default is 0.95).
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame retaining precursors with missingness below the selected threshold.
    module_output : str
        Updated log string.
    """
    timestamp = lambda: datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Validate allowed values for optimum_pr
    if precursor_thresh is not None:
        if not isinstance(precursor_thresh, (int, float)):
            module_output += f"[{timestamp()}] ###ERROR###: Invalid type for optimum_pr. Expected int or float, got {type(precursor_thresh).__name__}.\n"
            return data_df, module_output
        if not (0 <= precursor_thresh <= 100):
            module_output += f"[{timestamp()}] ###ERROR###: Invalid value for optimum_pr. Must be between 0 and 100. Got {precursor_thresh}.\n"
            return data_df, module_output

    # Validate allowed values for alpha_pr
    if not isinstance(precursor_alpha, (int, float)):
        module_output += f"[{timestamp()}] ###ERROR###: Invalid type for precursor alpha. Expected int or float, got {type(precursor_alpha).__name__}.\n"
        return data_df, module_output
    if not (0 <= precursor_alpha <= 1):
        module_output += f"[{timestamp()}] ###ERROR###: Invalid value for precursor alpha. Must be between 0 and 1. Got {precursor_alpha}.\n"
        return data_df, module_output
    
    # Calculate missingness and optimal threshold
    missingness, optimal_threshold = calculate_optimum_threshold(data_df.T, precursor_alpha, statistics='non-parametric')
    
    # Override with user-specified threshold if provided
    if precursor_thresh is not None:
        optimal_threshold = precursor_thresh / 100
        module_output += f"[{timestamp()}] Using hard thresholding at {precursor_thresh} % missingness for precursor filtering.\n"
    else:
        module_output += f"[{timestamp()}] Using automated thresholding {optimal_threshold:.2f} % missingness (alpha={sample_alpha}) for precursor filtering.\n"

    # Filter and return
    return data_df.loc[missingness < optimal_threshold], module_output

def filter_precursors_by_peptide_labels(data_df, module_output, remove_label=None):
    """
    Specify a e.g. SILAC label (part of the Precursor name), all precursors containing the label will be removed from further analysis
    author: BZ, SB
    
    Parameters
    ----------
    data_df : pd.DataFrame
        wide-format profiles, columns are samples
    module_output : str
        Log file to pick up run time comments and warnings
    remove_label : str
        name of the label to identify precursors to be removed, default = None

    Returns
    -------
    pd.DataFrame
        wide-format profiles, columns are samples, index with Precursor IDs.
    module_output : str
        Updated log string.
    """
    
    timestamp = lambda: datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    if not remove_label:
        module_output += f"[{timestamp()}] ###ERROR###: No precursors removed: remove_label is empty or None.\n"
        return data_df, module_output

    flagged_precursors = data_df.index.str.contains(remove_label)

    if flagged_precursors.any():
        module_output += f"[{timestamp()}] {flagged_precursors.sum()} precursors excluded because they contained the label: '{remove_label}'\n"
    else:
        module_output += f"[{timestamp()}] No precursors excluded: label '{remove_label}' not found in any index.\n"
    
    return data_df.drop(flagged_precursors, axis=0), module_output

def filter_proteins_by_peptide_count(data_df, id_mapping_df, module_output, min_peptides=1, protein_col="Protein.Group", precursor_col="Precursor.Id"):
    """
    Remove proteins identified by less precursors than specified with min_peptides.
    author: BZ, SB
    
    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    id_mapping_df: pd.DataFrame
        dataframe that stores mapping of precursors and proteins. 
    module_output : str
        Log file to pick up run time comments and warnings
    min_peptides : int
        minimum number of precursors required per protein, proteins with less precursors will be excluded. The default is 1. 
    protein_col : str
        Name of the ID column containing the protein or gene ID. The default is "Protein.Group".
    precursor_col : str
        Name of the ID column in id_mapping_df, containing the precursor ID. The default is "Precursor.Id".

    Returns
    -------
    data_df : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    module_output : str
        Updated log string.
    """
    # Count precursors / protein
    protein_counts = id_mapping_df[protein_col].value_counts()

    # Keep proteins with at least n >= min_peptdies
    proteins_to_keep = protein_counts[protein_counts >= min_peptides].index

    # Get precursors mapped to filtered proteins
    precursors_to_keep = id_mapping_df[id_mapping_df[protein_col].isin(proteins_to_keep)][precursor_col].drop_duplicates()

    # Filter on proteome level
    filtered_df = data_df.loc[data_df.index.isin(precursors_to_keep)]

    # Update logfile with dropped precursor count
    timestamp = lambda: datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    original_precursor_count = data_df.shape[0]
    filtered_precursor_count = filtered_df.shape[0]
    precursors_dropped = original_precursor_count - filtered_precursor_count
    module_output += (
        f"[{timestamp()}] " 
        f"Filtered out {int(precursors_dropped)} precursors "
        f"(kept {int(filtered_precursor_count)} out of {int(original_precursor_count)}) "
        f"based on minimum of {min_peptides} peptides per protein.\n"
    )
    return filtered_df, module_output

##############################################################################
#  Other Functions (Automated Thresholding, Label Replacement, Mask Updater)
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

def replace_labelled_precursors(data_df, module_output, target_label=None, replace_label=None): 
    """
    Specify a e.g. SILAC label (part of the Precursor name), all precursors containing the label will be replaced with another label
    author: BZ, SB
    
    Parameters
    ----------
    data_df : pd.DataFrame
        wide-format profiles, columns are samples
    target_label : str
        name of the label to identify precursors to be replaced. Default is None.
    replace_label : str
        Substring to replace the target_label with. Default is None.
    module_output : str
        logfile to update

    Returns
    -------
    pd.DataFrame
        wide-format profiles, columns are samples, index with Precursor IDs.
    module_output : str
        Updated log string.
    """
    timestamp = lambda: datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    if not target_label:
        module_output += f"[{timestamp()}]: ###ERROR###: No replacement performed: target_label is empty or None.\n"
        return data_df, module_output
    
    # Find matching precursors
    flagged_precursors = data_df.index.str.contains(target_label)

    # Perform replacement
    if flagged_precursors.any():
        # Replace in index
        new_index = data_df.index.to_series().replace(to_replace=target_label, value=replace_label, regex=True)
        # Update index
        data_df.index = new_index
        # Update log
        module_output += f"[{timestamp()}]: {flagged_precursors.sum()} precursors replaced: '{target_label}' -> '{replace_label}'\n"
    else:
        module_output += f"[{timestamp()}]: No replacement performed: {target_label} precursors not found\n"

    return data_df, module_output

def update_mask_df(data_df, data_processed_df, axis, method,
                   mask_df, mask_description_df, id_mapping_df, 
                   metadata_df, index_column):
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

    Parameters
    ----------
    data_df : pd.DataFrame
        The original dataset before filtering or processing.

    data_processed_df : pd.DataFrame
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

    id_mapping_df : pd.DataFrame
        Mapping DataFrame containing at least `Precursor.Id` column, 
        updated when `axis='precursors'`.

    metadata_df : pd.DataFrame
        Metadata associated with samples or precursors, updated to 
        reflect the new masking.

    index_column : str
        Column name in `metadata_df` to set as index for updates.

    Returns
    -------
    tuple of pd.DataFrame
        Updated `(mask_df, mask_description_df, metadata_df, id_mapping_df)`.

    Notes
    -----
    - When `axis='samples'`, masking is applied to columns present in `data_df` 
      but missing from `data_processed_df`.
    - When `axis='precursors'`, masking is applied to rows present in `data_df` 
      but missing from `data_processed_df`.
    - Both `mask_description_df` and `metadata_df` are reindexed before returning.
    """
        
    # Get the second largest value in mask_description_df safely
    sorted_values = mask_description_df['value'].sort_values()
    if len(sorted_values) >= 2:
        second_max_mapping = sorted_values.iloc[-2]
    else:
        second_max_mapping = sorted_values.iloc[-1]

    new_value = second_max_mapping + 1

    # Ensure method exists in mask_description_df or create it
    if method in mask_description_df.index:
        mask_description_df.loc[method, 'value'] = new_value
    else:
        mask_description_df.loc[method] = {'value': new_value}

    # Set index column for metadata_df
    metadata_df.set_index(index_column, inplace=True)

    # Update based on axis
    if axis == 'samples':
        flags = data_df.columns.difference(data_processed_df.columns)
        if not flags.empty:
            mask_df.loc[:, flags] = new_value
            metadata_df.loc[flags, method] = new_value

    elif axis == 'precursors':
        flags = data_df.index.difference(data_processed_df.index)
        if not flags.empty:
            mask_df.loc[flags, :] = new_value
            id_mapping_df.loc[id_mapping_df['Precursor.Id'].isin(flags), method] = new_value

    # Reset metadata_df and mask_description_df index for return
    metadata_df.reset_index(inplace=True)

    return mask_df, mask_description_df, metadata_df, id_mapping_df

##############################################################################
# Replicate Functions
##############################################################################

def filter_peptides_with_excessing_nans(peptide_profiles, file=None,profile_annotation_df=None,min_peptides=1,nan_threshold=1/3, mask=None, mask_description=None,values=None,ids=None,keep_track=True,diann_sample_column="File.Name"):
    """
    Remove precursors with too many missing values in given peptide_profiles.
    author: BZ

    # should work without multi-index 

    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    file : pd.DataFrame
        DIA-NN output report.tsv. The default is None.
    profile_annotation_df : pd.DataFrame
        Metadata mapping File.Name to shorter run names. The default is None.
    min_peptides : int
        minimum number of precursors required to quantify a protein. The default is 1.
    nan_threshold : float, optional
        Threshold for accepted number of missing values, all precursors with more than nan_threshold missing values in samples containing their respective protein are excluded from further analysis . The default is 1/3.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, generating the mask with the DIA-NN output file.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None.
    values : str
        values extracted from DIA-NN report.tsv, e.g. Precursor.Quantity. The default is None.
    ids : list, optional
        list of ids extracted from DIA-NN report.tsv. The default is ["Precursor.Id","Stripped.Sequence","Protein.Group","Genes"]. The default is None.
    keep_track : bool
        if set to True, keeps track of changes in filtering matrix (mask). The default is True.
    diann_sample_column = "File.Name" #name of the diann sample column to match with metadata, either "File.Name" or "Run".

    Returns
    -------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.

    """
    print("Excluding Precursors with excessive NaNs...")
    num_peptides=len(peptide_profiles)
    if global_variables.logfilename:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tRemoving Precursors with excessive missing values...\n")
        logfile.write("\t\t\t\tExcluding all values with more than "+str(nan_threshold)+" missing values in samples with presence of their protein.\n")
        logfile.close()
    if keep_track:    
        # if no mask, create mask with all zeros, except ID columns and sample group row
        # zero in mask is a value not excluded due to filtering, every other number describes reason of filtering, which is explained in mask_description
        if mask is None:
            if ids is None:
                print("Please provide Ids")
            mask,mask_description=module_helper.create_filtering_matrix(file, profile_annotation_df, min_peptides,values=values,ids=ids,diann_sample_column=diann_sample_column)
    if "sample group" in peptide_profiles.index:
        indices = peptide_profiles.index.drop("sample group")
    else:
        indices = peptide_profiles.index
    for pg in indices.get_level_values("Protein.Group").unique():
        temp_df = peptide_profiles[peptide_profiles.index.get_level_values("Protein.Group")==pg].replace(0, np.nan)
        indices_before_filtering = temp_df.index
        if len(temp_df)>1: # exclude fractions/samples without any precursor of this protein
            temp_df.dropna(thresh=int(np.ceil((1-nan_threshold)*temp_df.dropna(axis=1, how="all").shape[1])), inplace=True)  # keep
        else: #include all fractions/samples if protein is only represented by 1 precursor
            temp_df.dropna(thresh=int(np.ceil((1-nan_threshold)*temp_df.shape[1])), inplace=True)  # keep
        peptide_profiles.drop(indices_before_filtering.drop(temp_df.index), inplace=True)
        #remember in mask deleted peptides, but leave numbers marking not measured precursors
        if keep_track:
            #mask.loc[indices_before_filtering.drop(temp_df.index),peptide_profiles.columns]=mask.loc[indices_before_filtering.drop(temp_df.index),peptide_profiles.columns].replace(0,mask_description.loc["excessing_NaNs","value"])
            mask.loc[indices_before_filtering.drop(temp_df.index),peptide_profiles.columns]=mask.loc[indices_before_filtering.drop(temp_df.index),peptide_profiles.columns].replace(mask_description.loc["measured and kept"].values,mask_description.loc["excessing_NaNs","value"])
        
    print("Done")
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+str(np.round(100-len(peptide_profiles)/num_peptides*100,2))+"% of precursors excluded due to excessive number of missing values ("+str(num_peptides-len(peptide_profiles))+" out of "+str(num_peptides)+" precursors)\n")
        logfile.close()
    print(str(np.round(100-len(peptide_profiles)/num_peptides*100,2))+"% of precursors excluded due to excessive number of missing values")
    return([peptide_profiles, mask, mask_description])

def exclude_samples_with_too_many_NaNs(profiles,thresh=1/2,replace_zeros_with_nan=True):
    """
    Excludes samples (columns in profiles) with too many missing values.
    #thresh determines how many measured values have to be at least in data 

    # this works with data_df
    
    Parameters
    ----------
    profiles : pd.DataFrame
        rows equal protein groups or precursor IDs, columns correspond to samples and values in the dataframe to intensities.
    thresh : float, optional
        threshold on how many IDs relative to total amount of IDs are too many so that samples will be excluded. The default is 1/2.

    Returns
    -------
    profiles_trimmed : pd.DataFrame
        rows equal protein groups or precursor IDs, columns correspond to samples and values in the dataframe to intensities, samples with too many missing IDs excluded.
    profiles_deleted : pd.DataFrame
        rows equal protein groups or precursor IDs, columns correspond to samples and values in the dataframe to intensities, only samples with too many missing IDs.
        
    """
    num_samples_before=profiles.shape[1]
    if replace_zeros_with_nan:
        profiles_trimmed=profiles.copy().replace(0,np.nan)
    else:
        profiles_trimmed=profiles.copy()
    profiles_trimmed=profiles_trimmed.dropna(axis=1,thresh=int(thresh*len(profiles)))
    num_samples_after=profiles_trimmed.shape[1]
    print(str(num_samples_before-num_samples_after)+" samples excluded due to too many missing peptides (at least "+str(thresh*100)+" % of precursors measured over all samples required).")  #TODO print statements, ff
    if global_variables.logfilename:
        logfile=open(global_variables.logfilename,"a")
        logfile.write("\t\t\t\t"+str(num_samples_before-num_samples_after)+" samples excluded due to too many missing peptides (at least "+str(thresh*100)+" % of precursors measured over all samples required).\n")
        logfile.close()
    profiles_deleted=profiles[list(filter(lambda x: x not in profiles_trimmed.columns,profiles.columns))] 
    return([profiles_trimmed,profiles_deleted])

def exclude_samples_with_too_few_precursors(profiles,thresh=1/2,on_mad=True):
    """
    Excludes samples (columns in profiles) with too few IDs.

        # this works with data_df

    Parameters
    ----------
    profiles : pd.DataFrame
        rows equal protein groups or precursor IDs, columns correspond to samples and values in the dataframe to intensities.
    thresh : float, optional
        Threshold on how many IDs need to be at least in sample to be kept for further analysis, thresh * median(number of precursors over samples) determines how many precursors need to be at least in sample to be kept. The default is 1/2.
    on_mad : bool, optional
        if True, filters on median(number precursors) - thresh * MAD(number preecursors). The default is True.

    Returns
    -------
    profiles_trimmed : pd.DataFrame
        rows equal protein groups or precursor IDs, columns correspond to samples and values in the dataframe to intensities, samples with too few IDs excluded.
    profiles_deleted : pd.DataFrame
        rows equal protein groups or precursor IDs, columns correspond to samples and values in the dataframe to intensities, only samples with too few IDs.

    """
    num_samples_before=profiles.shape[1]
    num_precs=(profiles>0).sum()
    if on_mad:
        min_precs=np.nanmedian(num_precs)-(thresh*mad(num_precs))
    else:
        min_precs=thresh*np.nanmedian(num_precs)
    profiles_trimmed=profiles.copy().replace(0,np.nan).dropna(axis=1,thresh=int(min_precs))
    profiles_trimmed.replace(np.nan,0,inplace=True)
    num_samples_after=profiles_trimmed.shape[1]
    if on_mad: # print statements into log file as one string
        print(str(num_samples_before-num_samples_after)+" samples excluded due to too few measured peptides (at least median of all samples - "+str(thresh)+" times MAD precursors measured over all samples required).")
        print("Median number of precursors/proteins: "+str(np.median(num_precs)))
        print("Minimum number of precursors/proteins: "+str(np.median(num_precs)-(thresh*mad(num_precs))))
    else:
        print(str(num_samples_before-num_samples_after)+" samples excluded due to too few measured peptides (at least "+str(thresh*100)+" % of median precursor number measured over all samples required).")
    if global_variables.logfilename:
        logfile=open(global_variables.logfilename,"a")
        if on_mad:
            logfile.write("\t\t\t\t"+str(num_samples_before-num_samples_after)+" samples excluded due to too few measured peptides (at least median of all samples - "+str(thresh)+" times MAD precursors measured over all samples required).\n")
        else:
            logfile.write("\t\t\t\t"+str(num_samples_before-num_samples_after)+" samples excluded due to too few measured peptides (at least "+str(thresh*100)+" % of median precursor number measured over all samples required).\n")
        logfile.close()
    profiles_deleted=profiles[list(filter(lambda x: x not in profiles_trimmed.columns,profiles.columns))] 
    return([profiles_trimmed,profiles_deleted])

def exclude_samples_based_on_metadata(profiles,column="dropout/excl.",keep=np.nan):
    """
    exclude samples based on a value in a certain column in metadata, e.g. dropout/excl., keep is the value in that column to keep values

        # this works with data_df

    Parameters
    ----------
    profiles : pd.DataFrame
        rows equal protein groups or precursor IDs, columns correspond to samples and values in the dataframe to intensities.
    column : str, optional
        name of the multi column level in profiles with values defining which samples should be excluded/kept. The default is "dropout/excl.".
    keep : , optional
        defines which value in column defines which samples should be kept. The default is np.nan.

    Returns
    -------
    profiles_trimmed : pd.DataFrame
        rows equal protein groups or precursor IDs, columns correspond to samples and values in the dataframe to intensities, samples as defined excluded.
    profiles_deleted : pd.DataFrame
        rows equal protein groups or precursor IDs, columns correspond to samples and values in the dataframe to intensities, only excluded samples.

    """
    num_samples_before=profiles.shape[1]
    profiles_trimmed=profiles.xs(keep,axis=1,level=column,drop_level=False)
    num_samples_after=profiles_trimmed.shape[1]
    print(str(num_samples_before-num_samples_after)+" samples excluded based on metadata column '"+str(column)+"'.")
    if global_variables.logfilename:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(str(num_samples_before-num_samples_after)+" samples excluded based on metadata column '"+str(column)+"'.\n")
        logfile.close()
    profiles_deleted=profiles[list(filter(lambda x: x not in profiles_trimmed.columns,profiles.columns))] 
    return([profiles_trimmed,profiles_deleted])

##############################################################################
#  Main
##############################################################################

if __name__ == "__main__":
    
    #load 
    data_df, parameters_dict, id_mapping_df, seed, module_output, metadata_df, mask_df, mask_description_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = "Filtering wide", metadata=True, id_mapping=True, mask=True)

    #params
    filtering_wide_methods=parameters_dict['filtering_wide_methods']
    index_column=parameters_dict['index_column']
    key_column=parameters_dict['key_column']
    sample_thresh=parameters_dict['sample_thresh']
    sample_alpha=parameters_dict['sample_alpha']
    precursor_thresh=parameters_dict['precursor_thresh']
    precursor_alpha=parameters_dict['precursor_alpha']
    remove_label=parameters_dict['remove_label']
    min_peptides=parameters_dict['min_peptides']
    protein_col=parameters_dict['protein_col']
    precursor_col=parameters_dict['precursor_col']
    target_label=parameters_dict['target_label']
    replace_label=parameters_dict['replace_label']
    metadata_idx=index_column
    
    #helper functions
    timestamp = lambda: datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    def _unpack(result):
        """Allow called functions to return (df, log) or just df."""
        if isinstance(result, tuple) and len(result) == 2:
            return result[0], result[1]
        return result, ""

    def _append_log(logs_list, log_str):
        if isinstance(log_str, str) and log_str.strip():
            logs_list.extend([line for line in log_str.splitlines() if line.strip()])

    #run
    logs = []
    methods = filtering_wide_methods.replace(" ", "").split(";")
    data_processed_df = data_df.copy()

    if not methods:
        logs.append(f"[{timestamp()}] ###WARNING### No filtering methods provided. Skipping filtering pipeline.")
    else:
        for method in methods:
            logs.append(f"[{timestamp()}] Precursors x samples before {method}: {data_processed_df.shape}")

            # snapshot BEFORE this method
            pre_method_df = data_processed_df.copy()

            if method == "filter_samples_based_on_metadata":
                result = filter_samples_based_on_metadata(
                    data_processed_df, metadata_df, "", index_column, key_column
                )
                data_processed_df, log_str = _unpack(result)
                _append_log(logs, log_str)

                # compare before vs after THIS method (samples)
                mask_df, mask_description_df, metadata_df, id_mapping_df = update_mask_df(
                    pre_method_df, data_processed_df, 'samples', method,
                    mask_df, mask_description_df, id_mapping_df, metadata_df, metadata_idx
                )

            elif method == "filter_samples_by_detection_proportion":
                result = filter_samples_by_detection_proportion(
                    data_processed_df, "", sample_thresh, sample_alpha
                )
                data_processed_df, log_str = _unpack(result)
                _append_log(logs, log_str)

                # samples
                mask_df, mask_description_df, metadata_df, id_mapping_df = update_mask_df(
                    pre_method_df, data_processed_df, 'samples', method,
                    mask_df, mask_description_df, id_mapping_df, metadata_df, metadata_idx
                )

            elif method == "filter_precursors_by_detection_proportion":
                result = filter_precursors_by_detection_proportion(
                    data_processed_df, "", precursor_thresh, precursor_alpha
                )
                data_processed_df, log_str = _unpack(result)
                _append_log(logs, log_str)

                # precursors
                mask_df, mask_description_df, metadata_df, id_mapping_df = update_mask_df(
                    pre_method_df, data_processed_df, 'precursors', method,
                    mask_df, mask_description_df, id_mapping_df, metadata_df, metadata_idx
                )

            elif method == "filter_precursors_by_peptide_labels":
                result = filter_precursors_by_peptide_labels(
                    data_processed_df, "", remove_label
                )
                data_processed_df, log_str = _unpack(result)
                _append_log(logs, log_str)

                # precursors
                mask_df, mask_description_df, metadata_df, id_mapping_df = update_mask_df(
                    pre_method_df, data_processed_df, 'precursors', method,
                    mask_df, mask_description_df, id_mapping_df, metadata_df, metadata_idx
                )

            elif method == "filter_proteins_by_peptide_count":
                result = filter_proteins_by_peptide_count(
                    data_processed_df, id_mapping_df, "", min_peptides, protein_col, precursor_col
                )
                data_processed_df, log_str = _unpack(result)
                _append_log(logs, log_str)

                # precursors
                mask_df, mask_description_df, metadata_df, id_mapping_df = update_mask_df(
                    pre_method_df, data_processed_df, 'precursors', method,
                    mask_df, mask_description_df, id_mapping_df, metadata_df, metadata_idx
                )

            elif method == "replace_labelled_precursors":
                result = replace_labelled_precursors(
                    data_processed_df, "", target_label, replace_label
                )
                data_processed_df, log_str = _unpack(result)
                _append_log(logs, log_str)

                # precursors
                mask_df, mask_description_df, metadata_df, id_mapping_df = update_mask_df(
                    pre_method_df, data_processed_df, 'precursors', method,
                    mask_df, mask_description_df, id_mapping_df, metadata_df, metadata_idx
                )

            else:
                logs.append(f"[{timestamp()}] ###WARNING### Unknown filtering method: {method} — skipping.")

            logs.append(f"[{timestamp()}] Precursors x samples after  {method}: {data_processed_df.shape}")
   
    # finalize output

    module_output = "\n".join(logs)

    helper.close_module(data_processed_df, mask_df, mask_description_df,
                        log_file, module_output, out_dir_fp, outname, True,)