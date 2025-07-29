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
    - update_mask_df: Updates mask_df with a unique identifier based on the before / after a specific processing step.
"""

# %% Imports
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

    If the specified key_column is not found after indexing metadata, logs a ERRPR: COLUMN KEY ERROR, 
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
        return data_df

    meta_df_indexed = meta_df.set_index(index_column)

    if key_column not in meta_df_indexed.columns:
        module_output+=f"[{timestamp()}] ###ERROR###: COLUMN KEY ERROR - key_column '{key_column}' not found in metadata. Returning unfiltered data.\n"
        return data_df
    
    module_output+=f"[{timestamp()}] All keys found ('{index_column}', '{key_column}') in metadata. Filtering wide by metadata.\n"
    # Fill NAs in key column with FALSE
    sample_flags = meta_df_indexed.loc[:, key_column].fillna(False) 
    # Sample interesection between meta_df and data_df
    sample_flags = sample_flags[sample_flags.index.isin(data_df.columns)] 
    return data_df.loc[:, ~sample_flags], module_output

def filter_samples_by_detection_proportion(data_df, module_output, sample_thresh=None, alpha_sample=0.95):
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
    alpha_sample : float 
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
    if not isinstance(alpha_sample, (int, float)):
        module_output += f"[{timestamp()}] ###ERROR###: Invalid type for alpha_sample. Expected int or float, got {type(alpha_sample).__name__}.\n"
        return data_df, module_output
    if not (0 <= alpha_sample <= 1):
        module_output += f"[{timestamp()}] ###ERROR###: Invalid value for alpha_sample. Must be between 0 and 1. Got {alpha_sample}.\n"
        return data_df, module_output

    missingness, optimal_threshold = calculate_optimum_threshold(data_df, alpha_sample)

    if sample_thresh is not None:
        optimal_threshold = sample_thresh / 100 

    return data_df.loc[:, missingness < optimal_threshold], module_output

##############################################################################
# Precursor and Protein Thresholding Functions
##############################################################################

def filter_precursors_by_detection_proportion(data_df, module_output, precursor_thresh=None, alpha_precursor=0.95):
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
    alpha_precursor : float
        Alpha (significance) level for threshold calculation between 0 and 1(default is 0.95).
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame retaining precursors with missingness below the selected threshold.
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
    if not isinstance(alpha_precursor, (int, float)):
        module_output += f"[{timestamp()}] ###ERROR###: Invalid type for alpha_pr. Expected int or float, got {type(alpha_precursor).__name__}.\n"
        return data_df, module_output
    if not (0 <= alpha_precursor <= 1):
        module_output += f"[{timestamp()}] ###ERROR###: Invalid value for alpha_pr. Must be between 0 and 1. Got {alpha_precursor}.\n"
        return data_df, module_output
    
    # Calculate missingness and optimal threshold
    missingness, optimal_threshold = calculate_optimum_threshold(data_df.T, alpha_precursor, statistics='non-parametric')
    
    # Override with user-specified threshold if provided
    if precursor_thresh is not None:
        optimal_threshold = precursor_thresh / 100

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

    """

    # Find matching precursors
    flagged_precursors = data_df.index.str.contains(remove_label)

    # Update log
    if flagged_precursors.any():
        timestamp = lambda: datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        module_output += f"[{timestamp()}] {flagged_precursors.sum()} precursors excluded because they contained the label: {remove_label}\n"

    # Drop rows
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

def replace_labelled_precursors(data_df, target_label, replace_label, module_output): 
    """
    Specify a e.g. SILAC label (part of the Precursor name), all precursors containing the label will be replaced with another label
    author: BZ, SB
    
    Parameters
    ----------
    data_df : pd.DataFrame
        wide-format profiles, columns are samples
    target_label : str
        name of the label to identify precursors to be replaced.
    replace_label : str
        Substring to replace the target_label with.
    Module_output : str
        logfile to update

    Returns
    -------
    pd.DataFrame
        wide-format profiles, columns are samples, index with Precursor IDs.

    """
    # Find matching precursors
    flagged_precursors = data_df.index.str.contains(target_label)

    # Perform replacement
    if flagged_precursors.any():
        # Replace in index
        new_index = data_df.index.to_series().replace(to_replace=target_label, value=replace_label, regex=True)
        # Update index
        data_df.index = new_index
        # Update log
        module_output += f"{flagged_precursors.sum()} precursors replaced: '{target_label}' -> '{replace_label}'\n"

    return data_df, module_output

def udpate_mask_df(data_df, data_processed_df, axis, mask_df, mask_dicitionary):
    """"""
    return mask_df

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
# %% MAIN
##############################################################################

if __name__ == "__main__":
    
    data_df, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = "Filtering wide", metadata=True, id_mapping=False, mask=True, return_data_fp=False)

    # Load parameters (from JSON later)
    # parameters_dict = ...

    # Initialize the working DataFrame (always in wide format: rows = precursors, columns = samples)
    data_processed_df = data_df.copy()
    
    # Load methods
    filtering_wide_methods = "filter_samples_based_on_metadata; filter_samples_by_detection_proportion; filter_precursors_by_detection_proportion; " \
    "filter_precursors_by_peptide_labels; filter_proteins_by_peptide_count"
    methods = filtering_wide_methods.split(";")
    timestamp = lambda: datetime.now().strftime('%Y-%m-%d %H:%M:%S')
     
    for method in methods:
        if method == "filter_samples_based_on_metadata":
            module_output += (f"[{timestamp()}] Precursors x samples before {method}: {data_processed_df.shape}\n")
            data_processed_df, module_output = filter_samples_based_on_metadata(data_processed_df, metadata_df, module_output) 
            #mask update
            module_output += (f"[{timestamp()}] Precursors x samples of data after {method}: {data_processed_df.shape}\n")

        elif method == "filter_samples_by_detection_proportion": 
            module_output += (f"[{timestamp()}] Precursors x samples of data before {method}: {data_processed_df.shape}\n")
            data_processed_df, module_output = filter_samples_by_detection_proportion(data_processed_df, module_output)
            #mask update
            module_output += (f"[{timestamp()}] Precursors x samples of data after {method}: {data_processed_df.shape}\n")

        elif method == "filter_precursors_by_detection_proportion": 
            module_output += (f"[{timestamp()}] Precursors x samples of data before {method}: {data_processed_df.shape}\n") 
            data_processed_df, module_output = filter_precursors_by_detection_proportion(data_processed_df, module_output)
            #mask update
            module_output += (f"[{timestamp()}] Precursors x samples of data after {method}: {data_processed_df.shape}\n")

        elif method == "filter_precursors_by_peptide_labels": 
            module_output += (f"[{timestamp()}] Precursors x samples of data before {method}: {data_processed_df.shape}\n") 
            data_processed_df, module_output = filter_precursors_by_peptide_labels(data_processed_df, module_output)
            #mask update
            module_output += (f"[{timestamp()}] Precursors x samples of data after {method}: {data_processed_df.shape}\n") 

        elif method == "filter_proteins_by_peptide_count": 
            module_output += (f"[{timestamp()}] Precursors x samples of data before {method}: {data_processed_df.shape}\n") 
            data_processed_df, module_output = filter_proteins_by_peptide_count(data_processed_df, id_mapping_df, module_output)
            #mask mupdate
            module_output += (f"[{timestamp()}] Precursors x samples of data after {method}: {data_processed_df.shape}\n") 

        else:
            module_output += (f"[{timestamp()}] ###WARNING### Unknown filtering method: {method} — skipping.\n")

    # Save the final processed data
    #output_path = os.path.join(out_dir_fp, 'data_filtered_wide.tsv')
    #data_processed_df.to_csv(output_path, sep='\t')
    #print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Filtered data saved to: {output_path}") 
    helper.close_module(data_processed_df, masked_df, log_file, module_output, out_dir_fp, outname, mask=True)
