"""
Input: wide-format data
Output: Summarised wide-format data

    - print_args: decorator that prints arguments of a function (sanity checks)
    - calc_ref_intensity: calculates reference intensities
    - apply_protein_filter: filter proteins by their number of runs that they are detected in.
    - remove_rows_cols: remove rows and columns from a matrix.
    - get_peptide_ratios: calculate peptide ratios between runs
    - get_reference_intensities: calculate reference intensities for each protein in each run
    - estimate_protein_profile: estimate protein quantities in each run
    - weighted_harmonic_mean: calculates weighted harmonic mean
    - geometric_mean: calculates the geometric mean of an array x
    - summarize_peptides: summarizes peptide intensities to protein intensities with own approach using a pythagorean mean on representative precursors for a protein
    - transform_profiles_to_longformat_table: reshapes wide format profiles with IDs as rows, samples as columns and abundances as values to long-format with ID, sample and values as columns
    - MaxLFQ_summarization: summarize precursors to protein intensities using MaxLFQ
    - precursor_sum: sums precursor intensities to obtain protein quantities

"""

# TODO
# Adjust doc-strings

#general
import pandas as pd
import numpy as np
import time
import functools
from ctypes import cdll, c_double, POINTER, c_int, util, CDLL
from itertools import combinations
#from pathlib import Path
from inspect import getsourcefile
from os.path import abspath
import os
#import platform
from numpy.ctypeslib import ndpointer
from tqdm import tqdm


#own scripts
#from Python.Snippets import global_variables as global_variables #kept as placeholder for logging until main is being constructed
import module_helper as helper 



def print_args(func):
    """ Decorator that prints arguments of a function (sanity checks).
    
    Authors: Franziska Kistner, Vadim Demichev
    
    """
    @functools.wraps(func)
    def wrapper(*args):
        print(f"Function '{func.__name__}' called with args: {args[1:]}")
        result = func(*args)
        return result

    return wrapper


def calc_ref_intensity(x, quantity):
    #Authors: Franziska Kistner, Vadim Demichev
    return np.log10(x.nlargest(3, quantity, keep='first')[quantity].sum())


@print_args
def apply_protein_filter(long_format_precursor_df, accession, n_runs=None):
    """
    Filter proteins by their number of runs that they are detected in.
    
    Authors: Franziska Kistner, Vadim Demichev
    
    Parameters
    ----------
    accession: str,
        Column name that contains the protein accession
    long_format_precursor_df: pandas.Dataframe
        Df that contains all precursors with protein accession number and a run identifier

    n_runs: int,
        Minimal number of runs

    Returns
    -------
    filter runs: list,
        List of proteins that fulfill the filter criteria

    """
    # Group by proteins
    grp = long_format_precursor_df.groupby(accession)
    filter_runs = grp["Run"].nunique().loc[lambda x: x >= n_runs].reset_index()[accession]
    # Filter to proteins which are differentially expressed
    print("-- done.")

    return filter_runs


def remove_rows_cols(mat, indices):
    """
    Remove rows and columns from a matrix.
    
    Authors: Franziska Kistner, Vadim Demichev
    
    Parameters
    ----------
    mat: np.array(),
        Matrix to remove rows and columns from
    indices: list,
        List of indices to remove

    Returns
    -------
    no_nan_mat: np.array(),
        Matrix with removed rows and columns

    """
    # Get the number of rows and columns in the matrix
    n_rows, n_cols = mat.shape

    # Create a boolean mask that is True for all rows/columns except the indices to remove
    mask = np.ones((n_rows, n_cols), dtype=bool)
    mask[indices, :] = False
    mask[:, indices] = False

    # Use the mask to select the desired rows and columns
    no_nan_mat = mat[mask].reshape(n_rows - len(indices), n_cols - len(indices))

    return no_nan_mat


def get_peptide_ratios(wide_format, cc, minimum_ratios=1):
    """ Calculate peptide ratios between runs.
    This function (excluding documentation and slight modifications) is taken from
    https://github.com/MannLabs/alphapept/blob/master/nbs/08_quantification.ipynb.
    
    Authors: Franziska Kistner, Vadim Demichev
    
    Parameters
    ----------
    wide_format: np.array()
        Precursor in wide format represent a pivoted version of peptide quantities in all runs.
    cc: array-like
        Any column combinations that give indices of an upper triangular matrix of dimension n_runs x n_runs.
    minimum_ratios: int
        Minimum number of ratios

    Returns
    -------
    r: np.array()
        Median peptide ratios across runs (only uppertriangular)

    """
    n_runs = wide_format.shape[1]
    # Replace 0 with NaNs
    wide_format[wide_format == 0] = np.nan
    # Initialize results array
    # Dimensions: n_runs x n_runs
    r = np.empty((n_runs, n_runs))
    r[:] = np.nan

    # Compute pairwise run ratios i.e. compute the ratio of quantities in every run to every other run
    # For speed reasons, we compute the upper triangular of this matrix only.
    for indices in cc:
        run1_idx = indices[0]
        run2_idx = indices[1]
        # Calculate the ratio profile i.e.,
        # Calculate the ratio of every peptide pair in the current run to the second run.
        ratio = wide_format[:, run1_idx] / wide_format[:, run2_idx]
        # Calculate the number of non_nan ratios
        # You can set a minimum number of ratios (len(ratios)). Default is 1.
        non_nan = np.sum(~np.isnan(ratio))
        if non_nan >= minimum_ratios:
            # Take of median of the ratio profile that is considered to be a proxy for the protein ratio.
            # Check Cox et al. for more information.
            ratio_median = np.nanmedian(ratio)
        else:
            ratio_median = np.nan

        r[run1_idx, run2_idx] = ratio_median
    # Set the diagonal values to be 1.
    np.fill_diagonal(r, 1)

    return r


def get_reference_intensities(df, quantity, accession, run):
    """
    Calculate reference intensities for each protein in each run.
    
    Authors: Franziska Kistner, Vadim Demichev
    
    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe that contains the protein accession, run identifier and the ms2-based quantity
    quantity: str
        Column name that contains the ms2-based quantity
    accession:
        Column name that contains the protein accession
    run: str
        Column name that contains the run identifier


    Returns
    -------

    """
    # Sort the DataFrame by the quantity in descending order
    df_sorted = df.sort_values(quantity, ascending=False)

    # Use groupby and head to get the top 3 rows per (protein, run) group
    top3_per_group = df_sorted.groupby([accession, run]).head(3)

    # Group by protein_accession_column and run_id_column again and calculate the sum of ms2_based_quantity
    sum_top3_per_group = top3_per_group.groupby([accession, run])[quantity].sum()

    # Take the log10 of the sums
    log10_sum_top3 = np.log10(sum_top3_per_group)

    # Create the reference DataFrame using pivot_table
    reference = log10_sum_top3.reset_index().pivot_table(index=accession, columns=run, values=quantity).rename_axis(
        None, axis=1)

    return reference


#@print_args
def estimate_protein_profile(precursor_df, expected_n_runs, run, accession, sequence, ms2_quantity, ms1_quantity=None, metric=None,
                             extend=False, threshold=None, libname="maxlfq.so"):
    """ Estimate protein quantities in each run. This implementation uses the MaxLFQ algorithm (Cox et al., 2014)
    including a C++ (maxlfq.cpp, has to be in the same directory) program wrapped by ctypes.
    Make sure to create the required shared cpp library before using this script with

    On Linux    
    $ g++ -c -fPIC maxlfq.cpp -O3 -o maxlfq.o
    $ g++ -shared -Wl,-soname,maxlfq -o maxlfq.so -fPIC maxlfq.cpp
    
    On Windows
    $ g++ -c -fPIC maxlfq.cpp -O3 -o maxlfq.o
    $ g++ -shared -o maxlfq.dll -Wl,--out-implib,maxlfq.a -Wl,--export-all-symbols -Wl,--enable-auto-image-base maxlfq.o 
    
    Authors: Franziska Kistner, Vadim Demichev
    
    Parameters
    ----------
    precursor_df: pandas.DataFrame
        Df that contains all precursors with protein accession number, run identifier and quantity
    expected_n_runs: int
        Number of runs
    run: str
        Name of the column that contains a run identifier
    accession: str
        Name of the column that contains the protein accession number
    sequence: str
        Name of the column in the dataframe that contains the precursor sequence
    ms2_quantity: str
        Name of the column that contains the (precursor) quantities
    ms1_quantity: str
        Name of the column that contains a MS1-based quantity
    metric: str,
        Name of the column that is used for filtering for high-quality MS1-based quantites
    extend: bool,
       Extend the MS2-based quantities by high-quality (threshold set) MS1-based quanities
    threshold: float,
        Metric threshold
    libname: str
        Name of the C++ shared library file

    Returns
    -------
    protein_estimate_across_runs:
        Protein quantities for all unique proteins across all runs. Proteins were sorted lexicographically.

    """
    
    # Ctypes initialization
    #path_libc = util.find_library("msvcrt")
    
    #if platform.system() == "Windows":    
    libname_trim = libname.split(".")[0]
    lib="\\".join(abspath(getsourcefile(lambda:0)).split("\\")[:-1]) + "\\" + libname
    lib_trim="\\".join(abspath(getsourcefile(lambda:0)).split("\\")[:-1]) + "\\" + libname_trim
    #lib=lib.replace("\\","/")
    try:
        cpp_library = cdll.LoadLibrary(lib_trim)
    except FileNotFoundError:
        # Ctypes initialization
        '''
        for cpp_lib in ["libgcc_s_seh-1","libstdc++-6","libgcc_s_dw2-1","libgccjit-0","libwinpthread-1"]:
            path_libcpp = util.find_library(cpp_lib)
            if path_libcpp is not None:
                CDLL(path_libcpp)
        '''
        path_to_deps = "\\".join(util.find_library("libgcc_s_seh-1").split("\\")[:-1])
        #os.environ['PATH'] = path_to_deps + os.pathsep + os.environ['PATH']
        os.add_dll_directory(path_to_deps)
        try:
            lib="\\".join(abspath(getsourcefile(lambda:0)).split("\\")[:-1]) + "\\" + libname
            lib_trim="\\".join(abspath(getsourcefile(lambda:0)).split("\\")[:-1]) + "\\" + libname_trim
            cpp_library = cdll.LoadLibrary(lib)
        except FileNotFoundError:
            try:
                lib="Python/Preprocessing/"+libname
                #util.find_library(libname_trim)
                cpp_library = cdll.LoadLibrary(lib)
            except FileNotFoundError:
                lib = util.find_library("\\".join(abspath(getsourcefile(lambda:0)).split("\\")[:-1]) + "\\" + libname_trim)
                cpp_library = cdll.LoadLibrary(lib)
            
    # Specify argument dtypes
    cpp_library.maxlfq.argtypes = [c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int]
    
    
    # Read dataframe
    # Required columns:
    # Run identifier, Protein accession number, precursor sequence, (normalised) peptide quantity
    if extend:
        df = precursor_df[[run, accession, sequence, ms2_quantity, ms1_quantity, metric]]
    else:
        df = precursor_df[[run, accession, sequence, ms2_quantity]]
    sorted_df = df.sort_values(by=[accession, run]).reindex()

    # Find all unique proteins
    # We will loop over these proteins to estimate the protein quantities per run
    unique_proteins = sorted_df[accession].unique()

    # Compute reference intensities
    # Sum up the quantity of the 3 most intense peptides for each run per protein
    print("Calculating reference intensities...")
    reference = get_reference_intensities(df, ms2_quantity, accession, run)
    print("--done.")

    # Compute occurrence matrix i.e.,
    # wide format of precursor intensities per run
    occurrence_mat = sorted_df.pivot(index=run, columns=sequence, values=ms2_quantity)

    if extend:
        print("Extend peptide quantities by high-quality MS1 quantities.")
        print(f"Threshold: {threshold}")
        print("-- done.")
        # Filter to high-quality (metric > thr) precursor
        high_quality_precursor = sorted_df[sorted_df[metric] > threshold]
        high_quality_precursor.loc[:, sequence] = high_quality_precursor[sequence] + str("_MS1")
        # Extend MS2 quantities with MS1-based quantities
        extension = high_quality_precursor.pivot(index=run, columns=sequence, values=ms1_quantity)
        # Update long format dataframe and occurrence matrix
        sorted_df = pd.concat([sorted_df, high_quality_precursor], ignore_index=True, axis=0)
        occurrence_mat = pd.concat([occurrence_mat, extension], axis=1)
    else:
        print("Peptide ratios are calculated based on MS2 quantities only.")

    # Initialize column combinations
    cc = [_ for _ in combinations(range(expected_n_runs), 2)]
    # Save peptides per protein in a map
    # The protein to peptide map is used to subset the dataframe
    prot2pep_map = sorted_df.groupby(accession)[sequence].unique().to_dict()

    # Results array initialization
    # protein_profile stores the protein quantities for all available runs
    # row: proteins column: runs
    protein_profile = np.full((len(unique_proteins), expected_n_runs), np.nan)

    print("MaxLFQ-ing...")
    # For each unique protein p
    for p, protein in enumerate(tqdm(unique_proteins)):
        # Filter to peptides of current protein p using the map that was created before
        curr_om = occurrence_mat[prot2pep_map[protein]]
        # Subset ref intensities to current protein
        curr_ref_intensities = np.array(reference.loc[protein].values, dtype=c_double)
        # Count number of non-nan values i.e., the number of runs that the protein was detected in
        n_runs = np.count_nonzero(~np.isnan(curr_ref_intensities))
        # Compute the peptide ratios between each run
        ratios = get_peptide_ratios(curr_om.values.T, cc)
        # C++ cannot handle NaNs, so we need to convert them to 0
        # C++ needs floats in log10
        ratios_ = np.nan_to_num(np.log10(ratios)).astype(c_double)

        # Remove nan values from reference intensities and ratios
        non_nan_idcs = np.where(~np.isnan(curr_ref_intensities))[0]

        if np.any(np.isnan(curr_ref_intensities)):
            nan_idcs = np.where(np.isnan(curr_ref_intensities))[0]
            curr_ref_intensities = curr_ref_intensities[non_nan_idcs]
            # Remove rows and columns from ratios matrix
            ratios_ = remove_rows_cols(ratios_, nan_idcs)

        # Specify results dtypes
        cpp_library.maxlfq.restype = ndpointer(dtype=c_double, shape=(n_runs,))

        # Initialize weights matrix
        # If you do not want your peptide ratios to be weighted, simply pass a matrix of 1 where ratios != 0
        # Set weights to 0 where ratios is 0
        weights_ = np.ones_like(ratios_, dtype=c_double)
        weights_[ratios_ == 0] = 0

        # Ctypes wrapper that passes reference intensities, peptide ratios and weights
        max_lfq_protein_quantity = cpp_library.maxlfq(int(0),
                                                      curr_ref_intensities.ctypes.data_as(POINTER(c_double)),
                                                      ratios_.ctypes.data_as(POINTER(c_double)),
                                                      weights_.ctypes.data_as(POINTER(c_double)),
                                                      int(n_runs))
        # Save the estimated protein quantity at position p in the results array
        # Note: the maxlfq script returns log10 quantities
        protein_profile[p, non_nan_idcs] = max_lfq_protein_quantity

    print("-- done.")
    
    #create dataframe 
    protein_lfq_df = pd.DataFrame(protein_profile, index=unique_proteins, columns=occurrence_mat.index)
    protein_lfq_df.index.name = accession
    
    return protein_lfq_df

####################
####################

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
    return(np.product(x)**(1/len(x)))
    


def topN_summarization(data_df, id_mapping_df, filter_mask, module_output, top_n_peptides=3,weights=None,pooling="mean",min_peptides=1, seq_col="Precursor.ID", protein_col = "Protein.Group"):
    """
    summarize peptide intensities to protein intensities
    top_n_peptides most abundant, highly correlated and least imputed precursors over samples are summarized to a protein intensity
    #ranks precursors on number of measured peptides (not imputed ones)
    #ranks precursors on abundance
    #ranks precursors on correlation
    all three ranks are normalized, their harmonic mean is calculated and the average as defined by pooling over top_n_peptides precursors is calculated to retrieve the protein intensity.
    pooling is how to summarize selected precursors to protein group profiles, one of pythagorean means:  {"mean","hmean","gmean"}
    
    author: BZ
    
    Parameters
    ----------
    data_df : pd.DataFrame (old: peptide_profiles)
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    id_mapping_df: pd.DataFrame
        information on mapping from precursor to protein (group) level. rows denote individual precursors
    #####filter_mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.
    module_output: str
        log file to be filled.
    top_n_peptides : int, optional
        Number of precursors to summarize to protein intensities. The default is 3.
    weights : iterable, optional
        Weights for calculation of harmonic mean over the scores. The default is None, assigning equal weights.
    pooling : str, optional
        one of {"mean","hmean","gmean"} to summarize selected precursors to proteins. "hmean" equals harmonic mean, "gmean" equals geometric mean and "mean" the arithmetic mean. The default is "mean".
    min_peptides : int, optional
        require at least min_peptides precursors to calculate a protein intensity, else the protein intensity is set to 0. The default is 1.
    seq_col : str
        name of index level containing precursor Ids. The default is 'Precursor.Id'.
    protein_col : str
        name of index level containing protein Ids. The default is 'Protein.Group'.
        
        

    Returns
    -------
    pg_profiles : pd.DataFrame
        rows equal protein groups, columns correspond to samples and values in the dataframe to intensities.
    filter_mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.

    """
    
    pg_profiles=pd.DataFrame()
    
    if weights is None:
        weights=[1,1,1]
    data_df.replace(np.nan,0,inplace=True)
    print("Begin Data Summarization: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
        
    module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tSummarizing precursor quantities to Protein.Group intensities using own method...\n"
    module_output += "\t\t\t\t- summarize top "+str(top_n_peptides)+" precursors by their "+pooling+"\n"
    module_output += "\t\t\t\t- for the summarization precursors were selected by high abundance, correlation and low number of imputed values with following weights: "+", ".join([str(x) for x in weights])+"\n"

    indices=data_df.index
    #filter_mask = pd.DataFrame(columns=["Precursor summarized for Protein"],index=indices, values=0)    
        
    # Check if columns are in id_mapping_df
    if protein_col not in id_mapping_df:
        module_output += "###ERROR### Group column for summariaztion not found"
        return data_df, id_mapping_df, module_output
    
    if seq_col not in id_mapping_df:
        module_output += "###ERROR### Sequence column for summariaztion not found"
        return data_df, id_mapping_df, module_output
    
    # add new column to id_mapping_df
    ##### ID Mapping as output for summarization
    ######## Add save function
    id_mapping_df["Summarisation"] = 0
        
    pgs_list = id_mapping_df[protein_col].unique()
    id_mapping_dict = {pg:list(id_mapping_df[id_mapping_df[protein_col]==pg][seq_col]) for pg in pgs_list}
    
    #for every protein group
    for pg in pgs_list:
        
        precursors_list = [precursor for precursor in id_mapping_dict[pg] if precursor in indices]
        
        if precursors_list:
            temp_df = data_df.loc[precursors_list]
            
            ##on weighted harmonic mean between normalized abundance, correlation and amount of measured samples
            hmean_df=pd.DataFrame(index=temp_df.index,columns=["norm. Abundance","norm. Correlation","norm. Imputation"])
            sums=temp_df.sum(axis=1)
            corrs=temp_df.T.corr()
            corrs_hmean=[weighted_harmonic_mean(corrs.loc[x]) for x in corrs.index]
            temp_mask=filter_mask.loc[temp_df.index].map(lambda x: int(str(float(x)).split(".")[-1]))
            measured_precs=(temp_mask!=0).sum(axis=1)
            
            hmean_df["norm. Abundance"]=sums/sums.max() #division by zero if peptide not abundant over all fractions => should be excluded in remove_peptides_with_excessing_nans
            hmean_df["norm. Correlation"]=corrs_hmean/np.max(corrs_hmean) #division by zero if all precursors correlate zero
            hmean_df["norm. Imputation"]=measured_precs/measured_precs.max() #division by zero if all values are imputed for precursor => should be excluded in remove_peptides_with_excessing_nans
            hmean_df["hmean"]=[weighted_harmonic_mean(hmean_df.loc[x],weights) for x in hmean_df.index]
            
            ##### Check what this coulmn is doing
            temp_df.loc[hmean_df.sort_values(by="hmean",ascending=False).index[:top_n_peptides]].sum()
        
            #set indices for filter_mask
            if pooling=="mean":
                mask_indices=temp_df.loc[hmean_df.sort_values(by="hmean",ascending=False).index[:top_n_peptides]].index
            else: 
                prec_df=temp_df.loc[hmean_df.sort_values(by="hmean",ascending=False).index[:top_n_peptides]]
                mask_indices=prec_df.index
               
            id_mapping_df_c = id_mapping_df.copy().set_index(seq_col)
            id_mapping_df_c.loc[mask_indices,"Summarisation"] = 1
            
            if len(temp_df)>=min_peptides:
                if pooling=="mean":
                    pg_profiles.loc[pg]=temp_df.loc[hmean_df.sort_values(by="hmean",ascending=False).index[:top_n_peptides]].mean()
                elif pooling=="hmean":
                    pg_profiles.loc[pg]=[weighted_harmonic_mean(prec_df[c].values.reshape(len(prec_df))) for c in prec_df.columns]
                elif pooling=="gmean":
                    pg_profiles.loc[pg]=[geometric_mean(prec_df[c].values.reshape(len(prec_df))) for c in prec_df.columns]
                #filter_mask.loc[mask_indices,data_df.columns]="s"
                
                
            #else: #protein excluded due to less than min_peptides precursors
            #    filter_mask.loc[mask_indices,data_df.columns]=filter_mask.loc[mask_indices,data_df.columns].replace(mask_description.loc["measured and kept"].values,mask_description.loc["Min. "+str(min_peptides)+" Peptides"])
            #    filter_mask.loc[mask_indices,data_df.columns]=filter_mask.loc[mask_indices,data_df.columns].replace(mask_description.loc["Imputed"].values,mask_description.loc["Min. "+str(min_peptides)+" Peptides"])
    module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tEnd data summarization.\n"
    pg_profiles.index.name="Protein.Group"
    #add gene names
    #trans_df=data_df.reset_index().set_index("Protein.Group")
    #trans_df=trans_df[~trans_df.index.duplicated(keep="first")] #remove duplicated PGs
    #pg_profiles["Genes"]=trans_df.loc[pg_profiles.index,"Genes"]
    #pg_profiles=pg_profiles.reset_index().set_index(["Genes","Protein.Group"])
    #pg_profiles.columns=pg_profiles.columns.get_level_values(0)
    #pg_profiles.columns.name=None
    #return(pg_profiles,filter_mask,mask_description,module_output)

    return(pg_profiles.copy(),id_mapping_df_c,module_output)
    
    
    


def transform_profiles_to_longformat_table(profiles, id_mapping_df, protein_col, seq_col):
    """
    Reshapes wide format profiles with IDs as rows, samples as columns and abundances as values to long-format with ID, sample and values as columns.
    
    author: BZ
    
    Parameters
    ----------
    profiles : pd.DataFrame
        rows: precursors, columns: samples

    Returns
    -------
    lf_df : pd.DataFrame
        columns: [ID, sample, value]

    """
    '''
    ap_df=pd.DataFrame(columns=["precursor","protein_group","filename","prec_intensity"])
    ap_df["precursor"]=profiles.index.get_level_values("Precursor.Id").tolist()*len(profiles.columns)
    ap_df["protein_group"]=profiles.index.get_level_values("Protein.Group").tolist()*len(profiles.columns)
    ap_df["filename"]=np.repeat(profiles.columns.get_level_values("sample").tolist(),len(profiles))
    ap_df["prec_intensity"]=profiles.melt()["value"]
    '''
    
    #profiles.reset_index(inplace=True)
    lf_df=profiles.melt(ignore_index=False)
    lf_df.reset_index(inplace=True)
    lf_df.columns = [seq_col,"sample","value"]
    id_mapping_df_c = id_mapping_df.copy().set_index(seq_col)
    ### Add protein_col
    lf_df[protein_col] = [id_mapping_df_c.loc[lf_df.loc[index,seq_col],protein_col] for index in lf_df.index]
    return(lf_df, id_mapping_df_c)
    



#props to Franziska Kistner and Vadim Demichev for the MaxLFQ implementation+wrapper
def MaxLFQ_summarization(data_df, id_mapping_df, module_output, seq_col="Precursor.Id", protein_col="Protein.Group"):
    """
    Summarize precursors to protein intensities using MaxLFQ (Cox et al., 2014), implemented by Franziska Kistner with C++
    
    author: BZ
    
    Parameters
    ----------
    peptide_profiles : pandas.DataFrame
        rows: precursors, columns: samples
    unlog : bool
        If True, reverse log-transformation to basis 2. Default is False.
    seq_col : str
        name of index level containing precursor Ids. The default is 'Precursor.Id'.
    protein_col : str
        name of index level containing protein Ids. The default is 'Protein.Group'.
    run_col : str
        name of column level containing sample/run Ids. The default is 'sample'.
    Returns
    -------
    pd.DataFrame

    """
    module_output += "\n"+time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tSummarizing Precursors to Protein using MaxLFQ...\n"
        
    #remember all multi-level ids on protein/gene level
    #ids=peptide_profiles.index.names
    #protein_id_types=set(ids) & set(["Protein.Group","Protein.Ids","Protein.Names","Genes"])
    #gene_prot_df=pd.DataFrame(columns=list(protein_id_types))
    #for id_type in protein_id_types:
    #    gene_prot_df[id_type]=peptide_profiles.index.get_level_values(id_type)
    #gene_prot_df.drop_duplicates(inplace=True)
    #initiate variables for MaxLFQ
    n_runs=len(data_df.columns)
    #transform to long format
    data_df=helper.unlog_transform_intensities(data_df)
    lf_df, id_mapping_df_c = transform_profiles_to_longformat_table(data_df.copy())
    lf_df.replace(0,1,inplace=True) 
    #quantity_col = "value"
    #caclulate maxlfq
    #log10_protein_profiles = estimate_protein_profile_from_dataframe(precursor_df=lf_df, run_id_column=run_col, protein_accession_column=protein_col, precursor_seq_column=seq_col, quantity_column=quantity_col, n_runs=n_runs)
    log10_protein_profiles = estimate_protein_profile(precursor_df=lf_df, run="sample", accession=protein_col, sequence=seq_col, ms2_quantity="value", expected_n_runs=n_runs)
    data_prot_df=10**log10_protein_profiles
    #add all remaining IDs
    #data_df.reset_index(inplace=True)
    #gene_prot_df.set_index(protein_col,inplace=True)
    #protein_id_types.remove(protein_col)
    #for id_type in protein_id_types:
    #    protein_profiles[id_type] = gene_prot_df.loc[protein_profiles[protein_col], id_type].tolist()
    #protein_profiles.set_index(list(set(ids) & set(["Protein.Group","Protein.Ids","Protein.Names","Genes"])),inplace=True)
    #set proteins whose all precursors had zero values back to zero
    
    data_df.index = [id_mapping_df_c.loc[index,protein_col] for index in data_df.index]
    
    peptide_sums = data_df.groupby(by="index",level=0).sum()
    
    # Reorder to remove temporary replacement
    indices_prot = [index for index in peptide_sums.index if index in data_prot_df.index]
    columns_prot = [column for column in peptide_sums.columns if column in data_prot_df.columns]
    peptide_sums = peptide_sums.loc[indices_prot,columns_prot]
    
    #peptide_sums.columns=peptide_sums.columns.get_level_values(run_col)
    data_prot_df[peptide_sums==0]=0
    #add groups to protein profiles 
    #protein_profiles=protein_profiles[peptide_profiles.columns.get_level_values(run_col)]
    #protein_profiles.columns=peptide_profiles.columns
    #protein_profiles.columns.names=peptide_profiles.columns.names
    
    module_output += "\n"+time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tSummarization done.\n"
    return(data_prot_df, module_output)

def precursor_sum(data_df, id_mapping_df, module_output, seq_col="Precursor.Id", protein_col="Protein.Group"):
    """
    sums precursor intensities to obtain protein quantities.

    author: BZ/OL
    
    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    protein_levels : TYPE, optional
        name of index level(s) containing protein Ids. The default is ["Protein.Group","Genes"].

    Returns
    -------
    pd.DataFrame.

    """
    id_mapping_df_c = id_mapping_df.copy().set_index(seq_col)
    data_df.index = [id_mapping_df_c.loc[index,protein_col] for index in data_df.index]
    
    data_df = data_df.groupby(by="index",level=0).sum()
    
    module_output += "\n"+time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tProteins summarized to Precursors as sum of Precursor intensities.\n"
    
    return(data_df, module_output)




if __name__ == "__main__":
    
    data_df, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = "Summarisation", metadata=True, id_mapping=True, mask=True)
    
    if parameters_dict["summarization_method"]=="MaxLFQ":
        data_processed_df, module_output = MaxLFQ_summarization(data_df, id_mapping_df, module_output, 
                                                                seq_col=parameters_dict["seq_col"], 
                                                                protein_col = parameters_dict["protein_col"])
    elif parameters_dict["summarization_method"] in ["mean","gmean","hmean"]:
        data_processed_df, id_mapping_df, module_output = topN_summarization(data_df, id_mapping_df, masked_df, module_output, 
                                                                             top_n_peptides=parameters_dict["top_n_peptides"],
                                                                             weights=parameters_dict["weights"],
                                                                             pooling=parameters_dict["pooling"],
                                                                             min_peptides=parameters_dict["min_peptides"], 
                                                                             seq_col=parameters_dict["seq_col"], 
                                                                             protein_col = parameters_dict["protein_col"])
    elif parameters_dict["summarization_method"]=="sum":
        data_processed_df, module_output = precursor_sum(data_df, id_mapping_df, module_output, 
                                                         seq_col=parameters_dict["seq_col"], 
                                                         protein_col = parameters_dict["protein_col"])
    else:
        data_processed_df = data_df
        module_output += "###WARNING###: Summarization method not found. Skipped module."
        
    if parameters_dict["summarization_method"] in ["mean","gmean","hmean"]:
        id_mapping_df.to_csv(os.path.join(out_dir_fp,"id_mapping.tsv"),sep="\t")
        
    helper.close_module(data_processed_df, masked_df, log_file, module_output, out_dir_fp,outname, mask=True)