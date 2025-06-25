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


#general
import numpy as np
from collections import Counter
import time

#own scripts
from Python.Snippets import global_variables as global_variables #left this as a place holder for logging, when incorporating logging into main, these parts will be removed.
from modules_python import module_helper

##############################################################################
# Sample thresholding based on statistical inference on detection proportion
# Author: SB
##############################################################################

def condition_on_missingness(data_in):
    """
    Calculates the percentage of missing values for each precursor (column) in the input DataFrame.
    Adds a 'missingness' attribute sorted by missing percentage.

    author: SB

    Parameters
    ----------
    data_in : pd.DataFrame
        Input DataFrame with samples as rows and precursors as columns.

    Returns
    -------
    missingness : pd.Series
        Series indexed by precursor names, indicating percentage of missing values.
    """
    missing = data_in.isna().mean(axis=0) * 100
    missingness = missing.sort_values()
    return missingness


def calculate_optimum_threshold(data_in, alpha = .95):
    """
    Calculates the optimal missingness threshold based on a parametric (normal distribution)
    statistical approach. Provides visualizations (KDE and pseudo CDF plots) to assist threshold 
    determination.

    author: SB

    Parameters
    ----------
    data_in : pd.DataFrame
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

    missingness = condition_on_missingness(data_in)

    # Parametric stats
    mu = np.mean(missingness)
    sigma = np.std(missingness, ddof=1)
    z = stats.norm.ppf(alpha)  # One-sided Z value
    upper_bound = mu + z * sigma

    # Plot empirical KDE with normal overlay
    plt.figure(figsize=(6,5))
    sns.kdeplot(missingness, fill=True, label='Empirical KDE')

    # Plot normal PDF with same mean and std
    x_vals = np.linspace(min(missingness), max(missingness), 1000)
    normal_pdf = stats.norm.pdf(x_vals, mu, sigma)
    plt.plot(x_vals, normal_pdf, linestyle='--', label='Normal PDF', color='black')
    # Add vertical line at upper bound
    plt.axvline(upper_bound, color='red', linestyle='--', label=f'{int(alpha*100)}% Upper Bound = {upper_bound:.2f}%')
    # Labels and legend
    plt.xlabel('Sample Missingness (%)')
    plt.title(f'Missingness PDF with Parametric CI (alpha = {alpha})')
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # Plot pseudo CDF using cumulative retention of samples (not %)
    thresholds = np.linspace(0, 100, 21)
    recovered_counts = (missingness.values[None, :] <= thresholds[:, None]).sum(axis=1)
    plt.figure(figsize = (6,5))
    plt.plot(thresholds, recovered_counts)
    plt.title(f'Pseudo CDF using sample retention (alpha = {alpha})')
    plt.xlabel('Sample Missingness Threshold %')
    plt.axvline(upper_bound, linestyle = '--', color = 'red')
    plt.show()

    return missingness, upper_bound

def detection_control(data_in, optimum = 'Auto'):
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

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame retaining precursors below missingness threshold.
    """
    
    missingness, optimal_threshold = calculate_optimum_threshold(data_in)

         ## Check this bit as dtype might not work as intended
         ## Parsing of alpha might be included (missing)
         
    if optimum != 'Auto':
        optimal_threshold = optimum 

    data_in = data_in.loc[:,missingness<optimal_threshold]

    return data_in


def remove_peptides_with_excessing_nans(peptide_profiles, file=None,profile_annotation_df=None,min_peptides=1,nan_threshold=1/3, mask=None, mask_description=None,values=None,ids=None,keep_track=True,diann_sample_column="File.Name"):
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

    
def remove_proteins_with_too_few_precursors(profiles,mask,mask_description,min_peptides,protein_col="Protein.Group",keep_track=True):
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


    

def remove_labelled_peptides(peptide_profiles,remove_label,index_level="Precursor.Id"):
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

