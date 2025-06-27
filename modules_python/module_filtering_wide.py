"""
Input: wide-format data
Output: Filtered wide-format data

functions:
     - remove_peptides_with_excessing_nans: removes precursors with too many missing values
     - remove_proteins_with_too_few_precursors: removes proteins identified by less precursors than specified
     - assign_missingness: determines whether no precursor of the proteins is measured and if so sets its value to zero, else np.nan for imputation
     - remove_labelled_peptides: remove all precursors containing a certain e.g. SILAC label (part of the Precursor name) from further analysis
     - replace_labelled_peptides: replaces all occurences of a certain string, e.g. "(SILAC-K-L)" over all IDs with a given substitute string, e.g. K.
    
"""

# Imports
import numpy as np
import pandas as pd
import json
import os
import time
import argparse
import scipy.stats as stats
import module_helper as helper

#from Python.Snippets import global_variables as global_variables #left this as a place holder for logging, when incorporating logging into main, these parts will be removed.

##############################################################################
# Sample Thresholding
##############################################################################

def filter_samples_based_on_metadata(data_df, meta_df, column='none'):
    """
    Filters rows from the input data DataFrame based on missing values in a specific metadata column.

    If the specified metadata column does not exist, returns the original data DataFrame unchanged.

    Author: SB

    Parameters
    ----------
    data_df : pd.DataFrame
        DataFrame containing the data, with samples as rows and features/precursors as columns.

    meta_df : pd.DataFrame
        DataFrame containing metadata, indexed the same way as data_df.

    column : str, optional
        Column name in meta_df to check for missing values. Default is 'none'.

    Returns
    -------
    pd.DataFrame
        Filtered version of data_df where corresponding rows in meta_df[column] are missing (NaN).
    """
    if column not in meta_df.columns:
        return data_df

    return data_df.loc[meta_df[column].isna()]

def condition_on_detection_proportion(data_df):
    """
    Calculates the percentage of missing values for each precursor (column) in the input DataFrame.
    Adds a 'missingness' attribute sorted by missing percentage.

    author: SB

    Parameters
    ----------
    input_df : pd.DataFrame
        Input DataFrame with samples as rows and precursors as columns.

    Returns
    -------
    missingness : pd.Series
        Series indexed by precursor names, indicating percentage of missing values.
    """
    missing = data_df.isna().mean(axis=0) * 100
    missingness = missing.sort_values()
    return missingness

def calculate_optimum_threshold(data_df, alpha = .95):
    """
    Calculates the optimal missingness threshold based on a parametric (normal distribution)
    statistical approach. 

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

    # Parametric stats
    mu = np.mean(missingness)
    sigma = np.std(missingness, ddof=1)
    z = stats.norm.ppf(alpha)  # One-sided Z value
    upper_bound = mu + z * sigma

    # Plot empirical KDE with normal overlay
    #plt.figure(figsize=(6,5))
    #sns.kdeplot(missingness, fill=True, label='Empirical KDE')

    # Plot normal PDF with same mean and std
    #x_vals = np.linspace(min(missingness), max(missingness), 1000)
    #normal_pdf = stats.norm.pdf(x_vals, mu, sigma)
    #plt.plot(x_vals, normal_pdf, linestyle='--', label='Normal PDF', color='black')
    # Add vertical line at upper bound
    #plt.axvline(upper_bound, color='red', linestyle='--', label=f'{int(alpha*100)}% Upper Bound = {upper_bound:.2f}%')
    # Labels and legend
    #plt.xlabel('Sample Missingness (%)')
    #plt.title(f'{data_in.projectname}\nMissingness PDF with Parametric CI (alpha = {alpha})')
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(os.path.join(data_in.outer_path, data_in.projectname + str(alpha) + '_recovery_PDF.pdf'))
    #plt.show()
    
    # Plot pseudo CDF using cumulative retention of samples (not %)
    #thresholds = np.linspace(0, 100, 21)
    #recovered_counts = (missingness.values[None, :] <= thresholds[:, None]).sum(axis=1)
    #plt.figure(figsize = (6,5))
    #plt.plot(thresholds, recovered_counts)
    #plt.title(f'{data_in.projectname}\nPseudo CDF using sample retention (alpha = {alpha})')
    #plt.xlabel('Sample Missingness Threshold %')
    #plt.axvline(upper_bound, linestyle = '--', color = 'red')
    #plt.savefig(os.path.join(data_in.outer_path, data_in.projectname + str(alpha) + '_recovery_pseudoCDF.pdf'))
    #plt.show()

    return missingness, upper_bound

def filter_samples_on_detection_proportion(data_df, optimum = -3.14, alpha = 0.95):
    """
    Integrated pipeline that filters precursors from the dataset based on a calculated or 
    user-defined optimal missingness threshold.

    author: SB

    Parameters
    ----------
    data_in : pd.DataFrame
        Input DataFrame.
    optimum : str or float, optional
        Missingness threshold to use; 'Auto' calculates using `calculate_optimum_threshold`.
    alpha : float, optional
        Confidence interval alpha level for threshold calculation (default 0.95).

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame retaining precursors below missingness threshold.
    """
    
    missingness, optimal_threshold = calculate_optimum_threshold(data_df, alpha)
         
    if optimum != -3.14:
        optimal_threshold = optimum/100 

    return data_df.loc[:,missingness<optimal_threshold]

def exclude_samples_with_too_many_NaNs(profiles,thresh=1/2,replace_zeros_with_nan=True):
    """
    Excludes samples (columns in profiles) with too many missing values.
    #thresh determines how many measured values have to be at least in data 
    
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
    if on_mad:
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
# Precursor Thresholding
##############################################################################

def filter_peptides_with_excessing_nans(peptide_profiles, file=None,profile_annotation_df=None,min_peptides=1,nan_threshold=1/3, mask=None, mask_description=None,values=None,ids=None,keep_track=True,diann_sample_column="File.Name"):
    """
    Remove precursors with too many missing values in given peptide_profiles.
    author: BZ

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
    
def filter_proteins_with_too_few_precursors(profiles,mask,mask_description,min_peptides,protein_col="Protein.Group",keep_track=True):
    """
    Remove proteins identified by less precursors than specified with min_peptides.
    author: BZ
    
    Parameters
    ----------
    profiles : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, generating the mask with the DIA-NN output file.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None.
    min_peptides : int
        minimum number of precursors required per protein, proteins with less precursors will be excluded.
    protein_col : str, optional
        Name of the ID-column containing the protein or gene ID. The default is "Protein.Group".
    keep_track : bool, optional
        if set to True, keeps track of changes in filtering matrix (mask). The default is True.

    Raises
    ------
    KeyError
        Raised if protein_col not in profiles.index.

    Returns
    -------
    profiles : pd.DataFrame
        rows equal precursor or protein ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, generating the mask with the DIA-NN output file.

    """
    if protein_col not in profiles.index.names:
        raise KeyError(protein_col +"not in specified ids")
    #count precursors per proteins
    protcounts=Counter(profiles.index.get_level_values(protein_col))
    #identify proteins with too few precursors (<min_peptides)
    prots_too_few_precs=list(dict(filter(lambda x: x[1]<min_peptides, protcounts.items())).keys())
    #exclude proteins with too few precursors
    profiles=profiles.drop(prots_too_few_precs,axis=0,level=protein_col)
    #remember them in mask
    if keep_track:
        slice_ids=tuple(prots_too_few_precs if x==protein_col else slice(None) for x in mask.index.names)
        mask.loc[slice_ids,profiles.columns]=mask.loc[slice_ids,profiles.columns].applymap(lambda x: eval("x+0."+str(mask_description.loc["Min. "+str(min_peptides)+" Peptides","value"])))
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+str(len(prots_too_few_precs))+" Proteins excluded due to less than "+str(min_peptides)+" precursors.\n")
        logfile.close()
    return([profiles,mask])

def assign_missingness(peptide_profiles,mask,mask_description,protein_col="Protein.Group",on_filter_mask=True,keep_track=True):
    """
    Determines whether no precursor of the proteins is measured and if so sets its value to zero.
    author: BZ
    
    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, generating the mask with the DIA-NN output file.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None.
    protein_col : string, optional
        name of the index level to groupby protein ids
    on_filter_mask : bool, optional
        flag specifying whether to treat previously filtered precursors as present. The default is True. 
    keep_track : bool, optional
        if True, remembers in mask filtered precursors. The default is True

    Returns
    -------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, generating the mask with the DIA-NN output file.

    """
    #remember previous sorting
    orig_indices=peptide_profiles.index
    #sort profiles by protein_col
    peptide_profiles.sort_index(axis=0,level=protein_col,inplace=True)
    #create mask defining missing at random or present (True) vs missing not at random (False) on protein-level
    if on_filter_mask: #filtered precursors count as presence of a protein
        missingness_mask=(mask>=1).groupby(level=protein_col).sum()>0
        missingness_mask=missingness_mask.loc[peptide_profiles.index.get_level_values(protein_col).unique(),peptide_profiles.columns]
    else: #only consider precursors present after filtering
        missingness_mask=peptide_profiles.groupby(level=protein_col).sum()>0
    #replace all missing values with -1
    peptide_profiles.replace(0,-1,inplace=True)
    #multiply values with mask, resulting in 0 for values missing not at random, -1 for values missing at random and >0 for measured values
    peptide_profiles=peptide_profiles.mul(missingness_mask,axis=0,level=protein_col)
    #replace -1 values with nans for imputation
    peptide_profiles.replace(-1,np.nan,inplace=True)
    #remember missingness in filter mask
    if keep_track:
        description_value=eval("0."+str(mask_description.loc["missing not at random","value"]))
        mask.loc[peptide_profiles.index,peptide_profiles.columns]=mask.loc[peptide_profiles.index,peptide_profiles.columns]+(peptide_profiles==0)*description_value
    peptide_profiles=peptide_profiles.loc[orig_indices]
    return([peptide_profiles,mask])

def filter_labelled_peptides(peptide_profiles,remove_label,index_level="Precursor.Id"):
    """
    Specify a e.g. SILAC label (part of the Precursor name), all precursors containing the label will be removed from further analysis
    author: BZ
    
    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        wide-format profiles, columns are samples, index is a pd.MultiIndex with Precursor IDs.
    remove_label : str
        name of the label to identify precursors to be removed.
    index_level : str, optional
        name of the profiles' index level containing the precursor IDs with the remove_label substring. The default is "Precursor.Id".

    Returns
    -------
    pd.DataFrame
        wide-format profiles, columns are samples, index is a pd.MultiIndex with Precursor IDs.

    """
    i=peptide_profiles.index.names.index(index_level)
    all_exclude_precs=set()
    for rlabel in remove_label:
        exclude_precs=set(filter(lambda x: rlabel in x[i], peptide_profiles.index))
        all_exclude_precs=all_exclude_precs | exclude_precs
    print(str(len(all_exclude_precs)) + " Precursors excluded, because they contained following label: " + "; ".join(remove_label))
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+str(len(all_exclude_precs)) + " Precursors excluded, because they contained following label: " + "; ".join(remove_label) + "\n")
        logfile.close()
    return(peptide_profiles.drop(list(all_exclude_precs),axis=0))

def replace_labelled_peptides(peptide_profiles,replace_label,to_be_replaced_with,index_level="Precursor.Id"):
    """
    Specify a e.g. SILAC label (part of the Precursor name), all label occurences will be substituted. E.g. From K(SILAC-K-L) to K.
    author: BZ
    
    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        wide-format profiles, columns are samples, index is a pd.MultiIndex with Precursor IDs.
    remove_label : str
        name of the label to identify precursors to be removed. E.g. "K\(SILAC-K-L\)"
    to_be_replaced_with : str
        name of the substitute label. E.g. "K"
    index_level : str, optional
        name of the profiles' index level containing the precursor IDs with the remove_label substring. The default is "Precursor.Id".

    Returns
    -------
    pd.DataFrame
        wide-format profiles, columns are samples, index is a pd.MultiIndex with Precursor IDs.

    """
   
    i=peptide_profiles.index.names.index(index_level)
    peptide_profiles.index=peptide_profiles.index.set_levels(peptide_profiles.index.levels[i].str.replace(replace_label,to_be_replaced_with,regex=True),level=i)
    #precursors = peptide_profiles.index.get_level_values(index_level)
    #precursors = [re.sub(replace_label,to_be_replaced_with,precursor) for precursor in precursors]
    #duplicates = [i for i in set(precursors) if precursors.count(i) > 1]
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+ " Precursor label: " + replace_label + " replaced with "+ to_be_replaced_with + " for all precursors\n")        
        logfile.close()
    return(peptide_profiles)

##############################################################################
# MAIN
##############################################################################

if __name__ == "__main__":
    
    # Get data from wrapper by parsing flags
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--datapath", required=True, 
        help="Filepath to datafile that will be processed") 
    parser.add_argument(
        "-md", "--metapath", required=True, 
        help="Filepath to metadatafile") 
    parser.add_argument(
        "-o", "--outpath", required=True, 
        help="Filepath for outputfile") 
    #parser.add_argument(
    #    "-p", "--parampath", required=True, 
    #    help="Filepath to parameter json") 
    #parser.add_argument(
    #    "-md", "--metapath", required=True, 
    #    help="Filepath to metadatafile") 
    #parser.add_argument(
    #    "-ma", "--maskpath", required=True, 
    #    help="Filepath to maskfile") 
    #parser.add_argument(
    #    "-l", "--log_file", required=True, 
    #    help="Filepath to logfile") 
    #parser.add_argument(
    #    "-c", "--counter", required=True, 
    #    help="Module counter") 
    args = parser.parse_args()
    data_df_fp = args.datapath
    meta_df_fp = args.metapath
    out_dir_fp = args.outpath
    #param_fp = args.parampath
    #mask_fp = args.maskpath
    #log_file = args.log_file
    #counter = args.counter
    
    ###### To be modified
    outname = "data_filtered_wide"
    
    # Load parameter_file
    parameters_dict = {
            'filtering_wide_method1': 'metadata',
            'filtering_wide_method2': 'sample_detection',
            'filtering_wide_method3': ''
            }
    with open(os.path.join(out_dir_fp, 'parameters.json'), 'w') as f:
        json.dump(parameters_dict, f, indent = 4)

    with open(os.path.join(out_dir_fp, 'parameters.json')) as f:
        parameters_dict = json.load(f)

    # Load data
    data_df = pd.read_csv(data_df_fp, sep = '\t', index_col = 0)
    meta_df = pd.read_csv(meta_df_fp, sep = '\t', index_col = 0)
    #if parameters_dict["normalization_method"] not in ["delayed"]:
    #   data_df = pd.read_csv(data_df_fp,sep="\t",index_col=0)
    #else:
    #    data_arr = np.load(data_df_fp)
        
    # Define seed
    #seed = parameters_dict["seed"] + counter
        
    # Load mask and log file
    """
    mask = load_mask() # load mask (in parent directory)
    """
    
    #module_output = time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Normalization\n"
    #module_output += "Parameters\n------------\n"
    #for key,value in parameters_dict.items():
    #    module_output += str(key)+"\t\t\t"+str(value)+"\n"
    #module_output += "------------\n"  
    
    """
    # Helper DataFrames -> Maybe automate upper bit?
    log = helper.logger_add_parameters(log_file, parameters_dict)
    """
    
    data_processed_df = data_df.copy()

    # Process data
    if parameters_dict['filtering_wide_method1'] == "metadata":
        #module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Performing quantile normalization.\n"
        data_processed_df = filter_samples_based_on_metadata(data_df, meta_df)
        
    if parameters_dict['filtering_wide_method2'] == "sample_detection":
        #module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Performing quantile ranked normalization.\n"
        data_processed_df = filter_samples_on_detection_proportion(data_df)

    #module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Normalization finished\.n"

    # Save data

    data_processed_df.to_csv(os.path.join(out_dir_fp,outname+".tsv"),sep="\t")

    
    

