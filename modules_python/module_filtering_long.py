# -*- coding: utf-8 -*-"
"""
Input: Data import from software e.g. DIA-NN
Output: Filtered wide-format data

         
         
functions:
    - read_input: reads DIA-NN output report.tsv or MSFragger output
    - select_filtering_columns: provides based on applied filters a pre-selection of column names to read from the DIA-NN output report.tsv to save RAM
    - filtering: apply filtering on DIA-NN output report.tsv
    - filtering_peptide_specific: some optional additional filters based on peptide properties (so far only based on precursor charge state) => not used
    - create_profiles: converts (filtered) DIA-NN output report.tsv to wide-format profiles with IDs as rows and samples as columns (as specified by metadata) => Can be optimised in runtime
    - load_profiles: loads wide-format profiles if existent, else reads DIA-NN output report.tsv, filters it and creates wide-format profiles

#TODO
- create a function which checks beforehand for unique matching of DIA-NN runs filename and metadata filenames, e.g. by loading much smaller diann stats report filename column
- add filtering from MSFragger output

"""
# %%

from tkinter import Tk
from tkinter.filedialog import askopenfilename
import pandas as pd
import numpy as np
import os.path
import time
import logging
import seaborn as sns
import matplotlib.pyplot as plt
from modules_python import module_helper


##############################################################################
# Retention based sample thresholding (long)
##############################################################################

def calculate_retention_stats(df, step = 20):
    """ 
    Calculates number of recovered samples across incremental threshold intervals between 
    0 - 100% sample completeness.
    Author: SB
    
    Parameters:
    -----------
    df : pandas.DataFrame
        A DataFrame with rows representing samples and columns representing 
        protein features (or other measurements).
    
    step : int
        number of intervals for thresholding bteween 0 - 100% sample completeness,
        default = 20 

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the retention statistics with:
        - Samples: Number of retained samples at each threshold.
        - Precursors: Number of retained protein features.
        - Missing Sample Percentage: Average percentage of missing values.
    """
    summary = []
    # Remove rows that are entierely NaN
    data = df.dropna(axis = 0, how = 'all')

    # Calculate Retention Thresholds
    for x in range(step + 1):  
        gradient = x * (100 / step)
        select_sample = data.isnull().mean(axis=1) * 100
        sample_filtered = data.loc[select_sample <= gradient, :]
        if sample_filtered.empty:
            sample_miss = 0.0
        else:
            sample_miss = sample_filtered.isnull().mean(axis=0).mean() * 100
        summary.append((sample_filtered.shape[0], sample_filtered.shape[1], sample_miss))

    # Return retention data
    retention_stats = pd.DataFrame(summary, columns=['Samples', 'Precursors', 'Missing Sample Percentage'])
    return retention_stats

def find_optimal_threshold(retention_stats, step = 20):
    """ 
    Identifies the threshold at which missing sample percentage stabilizes 
    using first and second derivatives.
    Author: SB

    Parameters
    ----------
    retention_stats : pd.DataFrame
        DataFrame containing sample retention statistics.

    step : int, optional
        Number of intervals for thresholding sample completeness (default: 20).
    
    Returns
    -------
    float or None
        Optimal completeness threshold or None if no stabilization point is found.
    """
    # Set threshold multiplier and calculate retention stats
    xmult = 100 / step
    # Calculate missing threshold (retention) gradient minima
    samples = retention_stats['Missing Sample Percentage']
    threshold = retention_stats.index
    first_derivative = np.gradient(samples, threshold)
    second_derivative = np.gradient(first_derivative, threshold)
    zero_idx = np.where(np.isclose(second_derivative, 0, atol = 1e-6))[0] #zero grad
    if zero_idx.size > 0:
        return zero_idx[0] * xmult  # Convert to percentage scale
    return None

def apply_retention_filter(df, step=20, logfilename=None):
    """ 
    Applies thresholding on samples based on retention completeness using first and second derivatives.
    Author: SB

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame where rows are samples and columns are protein features.

    step : int, optional
        Number of intervals for thresholding sample completeness (default: 20).
    
    logfilename : str or None, optional
        If provided, logs progress to a specified file (default: None, no logging).

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame based on sample completeness.

    """
    if df.empty:
        logging.warning("DataFrame is empty. Returning unfiltered data.")
        return df

    logging.info("Calculating retention statistics...")
    retention_stats = calculate_retention_stats(df, step)

    logging.info("Identifying optimal threshold...")
    optimal_threshold = find_optimal_threshold(retention_stats, step)

    if optimal_threshold is not None:
        logging.info(f"Optimal threshold found at {optimal_threshold:.2f}% sample completeness.")

        # Plot retention curve
        sns.lineplot(x=retention_stats.index * (100 / step), y=retention_stats['Missing Sample Percentage'])
        plt.axvline(x=optimal_threshold, linestyle="--", color="red", label=f"Threshold: {optimal_threshold:.2f}%")
        plt.ylabel('Missing Sample Percentage')
        plt.xlabel('Sample Completeness (%)')
        plt.title('Retention Threshold Detection')
        plt.legend()
        plt.show()

        # Apply filtering
        filtered_df = df[df.isnull().mean(axis=1) * 100 <= optimal_threshold]
        logging.info(f"Dropped {df.shape[0] - filtered_df.shape[0]} samples.")

        # Write to log file if specified
        if logfilename:
            with open(logfilename, "a") as logfile:
                logfile.write(f"Optimal threshold: {optimal_threshold:.2f}% sample completeness.\n")
                logfile.write(f"Dropped {df.shape[0] - filtered_df.shape[0]} samples.\n")

        return filtered_df
    logging.warning("No  minimum found. Returning unfiltered data.")
    return df

##############################################################################
# File handling
##############################################################################


def read_input(filename,file_columns=None,logfilename=None):
    """
    Reads DIA-NN output report.tsv or MSFragger output, if specified only provided columns.
    author: BZ
    
    Parameters
    ----------
    filename : str, optional
        string specifying DIA-NN output path and filename (report.tsv). The default is None, opening a GUI to select the file.
    file_columns : list, optional
        list specifying which columns to read from DIA-NN output. The default is None, reading all columns.
    logfilename : str
        name (path + filename) of the logfile. The default is None, disabling logging into logfile.
    Returns
    -------
    DIA-NN output pandas.DataFrame

    """
    
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tLoading file...\n")
        logfile.close()
        
    #filesize = os.path.getsize(filename)
    #ram=psutil.virtual_memory()
    #if filesize<ram.available:
    # read tsv file
    file = pd.read_csv(filename, usecols=file_columns, sep="\t") # ,dtype={"Q.Value":float,"GG.Q.Value":float,"PG.Q.Value":float,"Protein.Q.Value":float,"Lib.Q.Value":float})
    
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tFile loaded successfully.\n")
        logfile.close()

    return(file)





def select_filtering_columns(filters,ids=["Precursor.Id", "Stripped.Sequence", "Protein.Group", "Genes"],values_col="Precursor.Quantity",input_sample_col="File.Name"):
    """
    with given filters dictionary, ids and values columns this function creates a list with column names of interest for extraction from DIA-NN output report.tsv
    author: BZ
    
    Parameters
    ----------
    filters : dict
        dictionary containing the used filters and filtering parameters.
    ids : list
        list of ids extracted from DIA-NN report.tsv. The default is ["Precursor.Id", "Stripped.Sequence", "Protein.Group", "Genes"].
    values_col : str
        values extracted from DIA-NN report.tsv, e.g. Precursor.Quantity.
    input_sample_col : str, optional
        name of the sample names in the DIA-NN output, DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    
    Returns
    -------
    file_columns : list
        list of e.g. DIANN output report.tsv column names to read.

    """
    #Only read necessary columns to save RAM (these computers are weak)
    filtering_columns=list()
    if filters["filter_proteotypic"]:
        filtering_columns.append("Proteotypic")
    if filters["filter_by_qvalue"]:
        filtering_columns.append("Q.Value")
    if filters["filter_by_global_qvalue"]:
        filtering_columns.append("Global.Q.Value")
    if filters["filter_by_gg_qvalue"]:
        filtering_columns.append("GG.Q.Value")
    if filters["filter_by_pg_qvalue"]:
        filtering_columns.append("PG.Q.Value")
    if filters["filter_by_global_pg_qvalue"]:
        filtering_columns.append("Global.PG.Q.Value")
    if filters["filter_by_p_qvalue"]:
        filtering_columns.append("Protein.Q.Value")
    if filters["filter_by_lib_qvalue"]:
        filtering_columns.append("Lib.Q.Value")
    if filters["filter_by_lib_pg_qvalue"]:
        filtering_columns.append("Lib.PG.Q.Value")
    if filters["filter_by_PEP"]:
        filtering_columns.append("PEP")
    if filters["filter_by_MScorr"] or filters["filter_by_MScorr_evidence_combined"]:
        filtering_columns.append("Ms1.Profile.Corr")
    if filters["filter_by_evidence"] or filters["filter_by_MScorr_evidence_combined"]:
        filtering_columns.append("Mass.Evidence")
    file_columns = [input_sample_col] + ids + [values_col] + filtering_columns    
    return(file_columns)




##############################################################################
# FILTERING
##############################################################################

# filter DIA-NN output (except on number of peptides), by default all filters set to False (no filtering)
def filtering(file, filters, ids=["Precursor.Id","Stripped.Sequence","Protein.Group","Genes"], profile_annotation_df=None, mask=None, mask_description=None, keep_track=True, values="Precursor.Quantity", filename_col="filename", diann_sample_column="File.Name", logfilename=None):
    """
    Filters rows in DIA-NN output report.tsv on filters as specified in filters by True.
    Number of filtered precursors in every step are saved to logfile.
    
    author: BZ
    
    Parameters
    ----------
    file : pd.DataFrame
        DIA-NN output report.tsv.
    filters : dict
        dictionary containing the used filters and filtering parameters.
    ids : list
        list of ids extracted from DIA-NN report.tsv. The default is ["Precursor.Id","Stripped.Sequence","Protein.Group","Genes"].
    profile_annotation_df : pd.DataFrame
        Metadata mapping File.Name to shorter run names. The default is None. 
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, calling creating_filtering_matrix() to create this variable.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None, calling creating_filtering_matrix() to create this variable.
    keep_track : bool, optional
        Flag to activate/deactivate tracking of filtered values in mask. The default is True.
    values : str
        values extracted from DIA-NN report.tsv. The default is "Precursor.Quantity".
    filename_col : str
        name of the column in profile_annotation_df specifying filenames (including path) as in DIA-NN report.tsv. The default is "filename".
    diann_sample_column : str, optional
        name of the sample names in the DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    logfilename : str or None, optional
        if not None, writes progress to this string specifying the name of the file for logging. The default is None, disabling logging of progress. 
    Returns
    -------
    file : pd.DataFrame
        filtered DIA-NN output report.tsv.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.


    """
    len_file_before=len(file)
    #Exclude runs from file not present in metadata
    if profile_annotation_df is not None:
        intersect_files=file[diann_sample_column].isin(profile_annotation_df[filename_col])
        if intersect_files.sum()==0:
            intersect_files=file[diann_sample_column].isin(profile_annotation_df.index)
        if intersect_files.sum()==0:
            print("WARNING: No overlap between samples in metadata and DIA-NN report. please check wheter filename_col is defined correctly.")
            return
        file=file[intersect_files]
        if len_file_before>len(file):
            print("Not all files in DIA-NN output are found in metadata. Analysis performed only on files present in metadata.")
            if logfilename is not None:
                logfile=open(logfilename,"a")
                logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tNot all files in DIA-NN output are found in metadata. Analysis performed only on files present in metadata.\n")
                logfile.close()
        if keep_track:
            if mask is None or mask_description is None:
                mask,mask_description=module_helper.create_filtering_matrix(file,profile_annotation_df,min_peptides=filters["min_peptides"],values=values,ids=ids,filename_col=filename_col,diann_sample_column=diann_sample_column)
    # FILTERING
    # FILTER BY PROTEOTYPIC
    # exclude non-proteotypic
    print("Filtering...")
    filtering_tracker_df=pd.DataFrame(index={key:value for key,value in filters.items() if value==True}.keys(),columns=["number of excluded precursors","number of remaining precursors"])
    len_file_before=len(file)
    len_file_after=len(file)
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tBegin Filtering...\n")
        logfile.close()
    if filters["filter_proteotypic"]:
        if "Proteotypic" in file.columns:
            if keep_track:
                rows=list(file[file["Proteotypic"]!=1][ids].itertuples(index=False))
                cols=file[file["Proteotypic"]!=1][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["Proteotypic","value"]
            file = file[file["Proteotypic"] == 1]
            filtering_tracker_df.loc["filter_proteotypic","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_proteotypic","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'Proteotypic' not in DIA-NN output, filter skipped")
    # FILTER BY Q-VALUE
    # exclude by q-value threshold for gene group AND precursor
    if filters["filter_by_qvalue"]:
        if "Q.Value" in file.columns:
            # replace xx,xxE+xx wich are seen as strings
            #file["Q.Value"]=list(map(lambda x: eval(x),file["Q.Value"]))
            if keep_track:
                rows=list(file[file["Q.Value"]>filters["qvalue_thresh"]][ids].itertuples(index=False))
                cols=file[file["Q.Value"]>filters["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["Q.Value","value"]
            file = file[file["Q.Value"] <= filters["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_q_value","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_q_value","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'Q.Value' not in DIA-NN output, filter skipped")
    if filters["filter_by_global_qvalue"]:
        if "Global.Q.Value" in file.columns:
            if keep_track:
                rows=list(file[file["Global.Q.Value"]>filters["qvalue_thresh"]][ids].itertuples(index=False))
                cols=file[file["Global.Q.Value"]>filters["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["Global.Q.Value","value"]
            file = file[file["Global.Q.Value"] <= filters["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_global_qvalue","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_global_qvalue","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'Global.Q.Value' not in DIA-NN output, filter skipped")
    if filters["filter_by_gg_qvalue"]:
        if "GG.Q.Value" in file.columns:
            if keep_track:
                rows=list(file[file["GG.Q.Value"]>filters["qvalue_thresh"]][ids].itertuples(index=False))
                cols=file[file["GG.Q.Value"]>filters["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["GG.Q.Value","value"]
            file = file[file["GG.Q.Value"] <= filters["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_gg_qvalue","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_gg_qvalue","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'GG.Q.Value' not in DIA-NN output, filter skipped")
    if filters["filter_by_pg_qvalue"]:
        if "PG.Q.Value" in file.columns:
            if keep_track:
                rows=list(file[file["PG.Q.Value"]>filters["qvalue_thresh"]][ids].itertuples(index=False))
                cols=file[file["PG.Q.Value"]>filters["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["PG.Q.Value","value"]
            file = file[file["PG.Q.Value"] <= filters["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_pg_qvalue","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_pg_qvalue","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'PG.Q.Value' not in DIA-NN output, filter skipped")
    if filters["filter_by_p_qvalue"]:
        if "Protein.Q.Value" in file.columns:
            if keep_track:
                rows=list(file[file["Protein.Q.Value"]>filters["qvalue_thresh"]][ids].itertuples(index=False))
                cols=file[file["Protein.Q.Value"]>filters["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["Protein.Q.Value","value"]
            file = file[file["Protein.Q.Value"] <= filters["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_p_qvalue","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_p_qvalue","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'Protein.Q.Value' not in DIA-NN output, filter skipped")
    if filters["filter_by_global_pg_qvalue"]:
        if "Global.PG.Q.Value" in file.columns:
            if keep_track:
                rows=list(file[file["Global.PG.Q.Value"]>filters["qvalue_thresh"]][ids].itertuples(index=False))
                cols=file[file["Global.PG.Q.Value"]>filters["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["Global.PG.Q.Value","value"]
            file = file[file["Global.PG.Q.Value"] <= filters["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_global_pg_qvalue","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_global_pg_qvalue","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'Global.PG.Q.Value' not in DIA-NN output, filter skipped")
    if filters["filter_by_lib_qvalue"]:
        if "Lib.Q.Value" in file.columns:
            if keep_track:
                rows=list(file[file["Lib.Q.Value"]>filters["lib_qvalue_thresh"]][ids].itertuples(index=False))
                cols=file[file["Lib.Q.Value"]>filters["lib_qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["Lib.Q.Value","value"]
            file = file[file["Lib.Q.Value"] <= filters["lib_qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_lib_qvalue","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_lib_qvalue","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'Lib.Q.Value' not in DIA-NN output, filter skipped")
    if filters["filter_by_lib_pg_qvalue"]:
        if "Lib.PG.Q.Value" in file.columns:
            if keep_track:
                rows=list(file[file["Lib.PG.Q.Value"]>filters["lib_qvalue_thresh"]][ids].itertuples(index=False))
                cols=file[file["Lib.PG.Q.Value"]>filters["lib_qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["Lib.PG.Q.Value","value"]
            file = file[file["Lib.PG.Q.Value"] <= filters["lib_qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_lib_pg_qvalue","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_lib_pg_qvalue","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'Lib.PG.Q.Value' not in DIA-NN output, filter skipped")

    # FILTER BY PEP
    # exclude rows with posterior error probability (PEP) for precursor identification exceeding pep_thresh
    if filters["filter_by_PEP"]:  # NOT CHECKED
        if "PEP" in file.columns:
            if keep_track:
                rows=list(file[file["PEP"]>filters["pep_thresh"]][ids].itertuples(index=False))
                cols=file[file["PEP"]>filters["pep_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["PEP","value"]
            file = file[file["PEP"] <= filters["pep_thresh"]]
            filtering_tracker_df.loc["filter_by_PEP","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_PEP","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Column 'PEP' not in DIA-NN output, filter skipped")
    # FILTER ON MS-LEVEL
    # exclude rows with correlation between MS1 and MS2 chromatograms below MScorr_thresh
    if filters["filter_by_MScorr"]:
        if "Ms1.Profile.Corr" in file.columns:
            if keep_track:
                rows=list(file[file["Ms1.Profile.Corr"] < filters["MScorr_thresh"]][ids].itertuples(index=False))
                cols=file[file["Ms1.Profile.Corr"] < filters["MScorr_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["Ms1.Profile.Corr","value"]
            file = file[file["Ms1.Profile.Corr"] >= filters["MScorr_thresh"]]
            filtering_tracker_df.loc["filter_by_MScorr","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_MScorr","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Ms1.Profile.Corr not in DIA-NN output, filter skipped")
    # exclude rows with MS2 evidence below evidence_thresh
    if filters["filter_by_evidence"]:
        if "Mass.Evidence" in file.columns:
            if keep_track:
                rows=list(file[file["Mass.Evidence"] < filters["evidence_thresh"]][ids].itertuples(index=False))
                cols=file[file["Mass.Evidence"] < filters["evidence_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask.loc[rows[r],cols[r]]=mask_description.loc["Mass.Evidence","value"]
            file = file[file["Mass.Evidence"] >= filters["evidence_thresh"]]
            filtering_tracker_df.loc["filter_by_evidence","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_evidence","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Mass.Evidence not in DIA-NN output, filter skipped")
    if filters["filter_by_MScorr_evidence_combined"]:  # NOT CHECKED
        if "Ms1.Profile.Corr" in file.columns and "Mass.Evidence" in file.columns:
            if keep_track:
                print("WARNING: tracking filtering for filter_by_MScorr_evidene_combined not implemented yet")#IMPLEMENT IF NECCESSARY
            file = file[file["Ms1.Profile.Corr"]+0.5 *
                        file["Mass.Evidence"] >= filters["MScorr_evidence_combined_thresh"]]
            filtering_tracker_df.loc["filter_by_MScorr_evidence_combined","number of excluded precursors"]=len_file_after-len(file)
            filtering_tracker_df.loc["filter_by_MScorr_evidence_combined","number of remaining precursors"]=len(file)
            len_file_after=len(file)
        else:
            print("Ms1.Profile.Corr and/or Mass.Evidence not in DIA-NN output, filter skipped")
    #filter out precursors with zero value
    if filters["filter_by_zero_values"]:
        if keep_track:
            rows=list(file[file[values] == 0][ids].itertuples(index=False))
            cols=file[file[values] == 0][diann_sample_column].tolist()
            for r in range(len(rows)):
                mask.loc[rows[r],cols[r]]=mask_description.loc["Zero "+values,"value"]
        file=file[file[values] > 0]
        filtering_tracker_df.loc["filter_by_zero_values","number of excluded precursors"]=len_file_after-len(file)
        filtering_tracker_df.loc["filter_by_zero_values","number of remaining precursors"]=len(file)
        len_file_after=len(file)
    print("Filtering Completed")
    len_file_after=len(file)
    filtering_tracker_df.dropna(inplace=True)
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tFiltering completed. "+str(len_file_before-len_file_after)+"("+str(np.round(100-len_file_after/len_file_before*100,2))+"%) Peptides were excluded during filtering of DIA-NN output\n\n")
        logfile.write(filtering_tracker_df.to_string().replace("\n","\n\t"))
        logfile.write("\n")
        logfile.close()
    if file.index.tolist() == []:
        print("Filtered file is empty. Please check if columns on which filters are based are not empty (sometimes depending on library used e.g. Genes and therefore Q-Values on gene-level are missing), else adjust filtering parameters.")
    return([file,mask,mask_description])


#Construction site, improve if needed
# filtering filters already peptide and thus protein specific, this function adds filters solely on peptide features
def filtering_peptide_specific(file, filter_mask, mask_description, filter_by_precursorcharge=False, preccharge_min=2, preccharge_max=42):
    """
    
    Additional filters that only apply on precursors.
    In development. Placeholder for filtering based on physicochemical peptide properties derived from Precursor/Peptide sequence, such as GRAVY, aliphatic index, peptide length, m/z, ...
    
    author: BZ
    
    Parameters
    ----------
    file : pd.DataFrame
        (filtered) DIA-NN output report.tsv.
    filter_mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. 
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.
    filter_by_precursorcharge : bool, optional
        flag to activate/deactivate filtering by Precursor.Charge. The default is False.
    preccharge_min : ing, optional
        minimum precursor charge required for a precursor to be considered for further analysis. The default is 2.
    preccharge_max : int, optional
        maximum precursor charge required for a precursor to be considered for further analysis. The default is 42.

    Returns
    -------
    file : pd.DataFrame
       filtered DIA-NN output report.tsv.

    """
    # filter on Precursor.Charge, exclude Precursors with charge below preccharge_min and higher preccharge_max
    if filter_by_precursorcharge:
        file = file[file["Precursor.Charge"] >= preccharge_min]
        file = file[file["Precursor.Charge"] <= preccharge_max]
        #keep track in filter_mask
        #...
    # e.g. remove low MS2 signal peptides
    # e.g. one could filter by parameters obtained with PeptideCalculations.py, such as GRAVY, InstabilityIndex etc.
    return(file)




#error in "Assigning measured proteins [...]" ff during function call if ids in profile_df does not match ids in filter_mask
def create_profiles(file, ids, values, profile_annotation_df, analysis_folder, savefilename, filter_mask=None, mask_description=None, min_peptides=1, count_peptides=True, filename_col="filename", sample_col="name on x-axis", sample_group_col="sample group", diann_sample_column="File.Name",logfilename=None):
    """
    Creates wide-format intensity profiles
    
    author: BZ
    
    Parameters
    ----------
    file : pd.DataFrame
        (filtered) DIA-NN output report.tsv.
    ids : list
        list of ids extracted from DIA-NN report.tsv.
    values : str
        values extracted from DIA-NN report.tsv.
    profile_annotation_df : pd.DataFrame
        Metadata mapping diann_sample_column (File.Name/Run) to shorter run names.
    analysis_folder : str
        path of analysis and where metadata should be saved.
    savefilename : str
        name of pre-processed file.
    filter_mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, calling creating_filtering_matrix() to create this variable.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None, calling creating_filtering_matrix() to create this variable.
    min_peptides : int
        minimum number of precursors required to quantify a protein. The default is 1.
    count_peptides : bool, optional
        flag to activate/deactivate counting of peptides. The default is True.
    filename_col : str, optional
        column name with complete filenames in metadata. The default is "filename".
    sample_col : str, optional
        column name with sample names in metadata. The default is "name on x-axis".
    sample_group_col: str, optional
        column name with sample group names in metadata. The default is "sample group". 
    diann_sample_column : str, optional
        name of the sample names in the DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    logfilename : str or None, optional
        if not None, writes progress to this string specifying the name of the file for logging. The default is None, disabling logging of progress. 
        
    Returns
    -------
    profile_df : pd.DataFrame.
        wide format intensities, rows correspond to ids, columns to filenames and intensities to values.
    num_peptide_df : pd.DataFrame or None
        wide format dataframe as profile_df with number of precursors for every id in ids and filenames, num_peptide_df is None if count_peptides=False.
    filter_mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, calling creating_filtering_matrix() to create this variable.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None, calling creating_filtering_matrix() to create this variable.
    

    """
    #write to logfile
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tCreating profiles...\n")
        logfile.close()
    # load filter_mask if needed and not provided
    if count_peptides and filter_mask is None:
        filter_mask,mask_description=module_helper.create_filtering_matrix(file,profile_annotation_df,min_peptides=min_peptides,values=values,ids=ids,diann_sample_column=diann_sample_column)
    # delete path from diann_sample_column and set index in file (DIA-NN output)
    file["FileName_trimmed"] = list(map(lambda x: x.replace("/", "\\").split("\\")[-1], file[diann_sample_column].tolist()))
    # exclude filenames not present in profile_annotation_df
    len_file_before=len(file)
    #file = file[list(map(lambda x: file.loc[x, "FileName_trimmed"] in profile_annotation_df["filename"], file.index))]
    file=file[[file.loc[x,"FileName_trimmed"] in profile_annotation_df[filename_col] for x in file.index]]
    if len_file_before>len(file):
        print("Not all files in DIA-NN output are found in metadata. Analysis performed only on files present in metadata.")
        if logfilename is not None:
            logfile=open(logfilename,"a")
            logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tNot all files in DIA-NN output are found in metadata. Analysis performed only on files present in metadata.\n")
            logfile.close()
    #file = file[list(map(lambda x: file.loc[x, "FileName_trimmed"] in profile_annotation_df.index, file.index))]
    # trim profile_annotation_df to intersecting raw files with file (columns pressent in resulting profiles_df)
    profile_annotation_df = profile_annotation_df.loc[list(set(
        file["FileName_trimmed"]) & set(profile_annotation_df.index))]
    profile_annotation_df.sort_values(
        by=[sample_group_col,sample_col], inplace=True, axis=0)
    #there should be no duplicates
    file=file.drop_duplicates()
    if type(ids) == list:  # e.g. ids=["Protein.Group","Genes"]
        file.set_index(["FileName_trimmed"] + ids, inplace=True)
    else:
        file.set_index(["FileName_trimmed", ids], inplace=True)
    file.sort_index(inplace=True)
    # check if provided ids column has unique entries (obligatory)
    if len(file.index.unique()) < len(file.index):
        print("WARNING: Multiple times occuring IDs found (ok on protein level)")
    print("Creating profiles...")
    #profile_df = pd.DataFrame(index=file.index.droplevel(0).unique(), columns=profile_annotation_df["name on x-axis"])
    profile_df = pd.DataFrame(index=file.index.droplevel(
        0).unique(), columns=pd.MultiIndex.from_arrays([profile_annotation_df[sample_col].tolist(),profile_annotation_df[sample_group_col].tolist()], names=('sample', 'group')))
    #if precursor/peptide intensities are of interest
    if "Precursor" in ",".join(ids) or "Stripped" in ",".join(ids):
        count_peptides = False
        index_dict={index[0]:index[0:] for index in file.index.droplevel(0).unique()}
        profile_df=pd.pivot_table(file.reset_index(),index='Precursor.Id',columns="FileName_trimmed",values=values)
        profile_df.columns= profile_df.columns.map(lambda x: tuple(profile_annotation_df.loc[x,[sample_col, sample_group_col]].tolist()))
        profile_df.index=profile_df.index.map(lambda x: index_dict[x])
        profile_df.columns=profile_df.columns.rename(['sample', 'group'])
        profile_df.index = profile_df.index.rename(ids)
        #for index in file.index:
        #    profile_df.loc[index[1:], tuple(profile_annotation_df.loc[index[0],
        #                                                        [sample_col,sample_group_col]].tolist())] = file.loc[index, values]#[0]
    else:
        for index in file.index:
            '''
            ID=index[1]
            xvalue=index[0]
            if np.round(file.loc[(xvalue, ID), values][0],0)!=np.round(file.loc[(xvalue, ID), values].mean(),0): 
                print("oh-oh, different values over peptides encountered")
            '''
            profile_df.loc[index[1:],tuple(profile_annotation_df.loc[index[0],
                                                                [sample_col,sample_group_col]].tolist())] = file.loc[index, values][0]  # [0] because maxLFQ.Unique s same over all peptides
    profile_df.replace(np.nan, 0, inplace=True)
    if count_peptides:
        num_peptide_df = profile_df.copy()
        for index in file.index:
            num_peptide_df.loc[index[1:], tuple(profile_annotation_df.loc[index[0],
                                                                [sample_col,sample_group_col]].tolist())] = len(file.loc[index, values])  # len(file.loc[(index[0], index[1]), values])
        # apply filter based on number of peptides
        profile_df = profile_df*(num_peptide_df >= min_peptides)
        # keep track of filtered values in filtering matrix
        #Implement if needed, filtering matrix has indices on precursor level, but num_peptide_df may be on protein level only
        for col in filter_mask.columns:
            filter_mask.loc[num_peptide_df[num_peptide_df[col]>=min_peptides].index,col]=mask_description.loc["Min. "+str(min_peptides)+" Peptides"] ###
    else:
        num_peptide_df = None
        # num_peptide_df=pd.DataFrame(index=profile_df.index,columns=profile_df.columns)
    # exclude all with only zeros in profile
    profile_df.drop(profile_df[profile_df.sum(axis=1) == 0].index, inplace=True)
    # SAVE
    profile_df.index.set_names(ids, inplace=True)
    profile_df.to_csv(analysis_folder +
                      savefilename+".tsv", sep="\t")
    filter_mask.to_csv(analysis_folder+savefilename+"_filtermask.tsv",sep="\t") #saved as backup, will be overwritten later
    mask_description.to_csv(analysis_folder+savefilename+"_filtermask_description.tsv",sep="\t")
    if count_peptides:
        num_peptide_df.replace(np.nan, 0, inplace=True)
        num_peptide_df.index.set_names(ids, inplace=True)
        num_peptide_df.to_csv(analysis_folder +
                              savefilename+"_num_peptides.tsv", sep="\t")
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tProfiles created and saved\n")
        logfile.close()
    #add to filter_mask precursors that were measured as 1
    print("Assigning measured precursors in filtering matrix...")
    for col in profile_df.columns:
        filter_mask.loc[profile_df[profile_df[col]>0].index,col]=mask_description.loc["measured and kept","value"]
    filter_mask.to_csv(analysis_folder+savefilename+"_filtermask.tsv",sep="\t") #saved as backup, will be overwritten in peptide-centric approach
    mask_description.to_csv(analysis_folder+savefilename+"_filtermask_description.tsv",sep="\t")
    return([profile_df, num_peptide_df, filter_mask, mask_description])



def load_profiles(ids, values=None, profiles_filename=None, raw_filename=None, filters=dict(), profile_annotation_file=None, profile_order=None, raw_data_folder="raw_data", profile_annotation="", filename2run_pattern=None, filename2samplegroup_pattern="1", overwrite_filename_annotations=False, filename_col="filename", encoding="utf-8", diann_sample_column="File.Name"):
    """
    tries to load the wide-format pd.DataFrame with ids' intensities, and if not found calls create_profiles to create the wide-format file.
    
    author: BZ
    
    Parameters
    ----------
    ids : list
        list of ids extracted from DIA-NN report.tsv.
    values : str
        values extracted from DIA-NN report.tsv.
    profiles_filename : str, optional
        name of the wide-format profiles pd.DataFrame created by create_profiles(). The default is None, calling a GUI to select the file.
    raw_filename : str, optional
        name of the DIA-NN output report.tsv file. The default is None, calling a GUI to select the file.
    filters : dict
        dictionary containing the used filters and filtering parameters. The default is dict().
    profile_annotation_file : str
        name of metadata file. The default is None, getting the 
    profile_order : list, optional
        list of File.Names/Runs how they should be ordered in metadata and thus in wide format profiles. The default is None.
    raw_data_folder : str, optional
        path to *.quant files, if provided for all files metadata will be generated, even if no proteins were detected in file. The default is None.
    profile_annotation : str, optional
        profile_annotation is one of ["manual","IEF","SEC","numbers","from_file",""]. "manual" requests input via console, "Excel" saves filenames from folder "raw" to a tsv file for further manual annotation, "IEF" requests input pH range and calculates according pH distribution, "SEC" requests input size range and calculates according distribution, "numbers" assigns integers starting from 1, "from file" opens GUI to select TSV file with annotated names, any other string leaves the filenames. The default is "", corresponding to copying the File.Name.
        Used in previous versions, not necessary anymore
    filename2run_pattern : str
        A with eval evaluable string to get a substring as shorter sample name from diann_sample_column, x is the iterable over all entities. The default is None, using the File.Name
    filename2samplegroup_pattern : str
        A with eval evaluable string to get a substring as shorter sample group from diann_sample_column, x is the iterable over all entities. The default is "1", assigning group 1 to every sample.
    overwrite_filename_annotations : bool, optional
        flag to activate/deactivate overwriting of existent automatically generated metadata. The default is False.
    filename_col : str
        name of the column in profile_annotation_df specifying filenames (including path) as in DIA-NN report.tsv. The default is "filename".
    encoding : str, optional
        encoding for metadata file passed to pandas.read_csv. The default is "utf-8".
    diann_sample_column : str, optional
        name of the sample names in the DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    
    Returns
    -------
    profiles_file : pd.DataFrame.
        wide format intensities, rows correspond to ids, columns to filenames and intensities to values.
    profiles_file_numpeptides : pd.DataFrame or None
        wide format dataframe as profile_df with number of precursors for every id in ids and filenames, num_peptide_df is None if count_peptides=False.
    filter_mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, calling creating_filtering_matrix() to create this variable.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None, calling creating_filtering_matrix() to create this variable.


    """
    if profiles_filename is None:
        print("Please select filtered profiles from DIA-NN output:")
        root = Tk()
        root.withdraw()
        profiles_filename = askopenfilename(initialdir="../")
        root.destroy()
    analysis_folder = "/".join(profiles_filename.split("/")[:-1])+"/"
    savefilename = profiles_filename.split("/")[-1]
    try:
        profiles_file = pd.read_csv(profiles_filename, sep="\t", index_col=list(range(len(ids))),header=[0,1])
        #profiles_file.set_index(ids, inplace=True)
        if "Precursor" not in ",".join(ids) or "Stripped" not in ",".join(ids):
            profiles_file_numpeptides = pd.read_csv(
                profiles_filename[:-4]+"_num_peptides.tsv", sep="\t",index_col=list(range(len(ids))),header=[0,1])
        else:
            profiles_file_numpeptides = None
        filter_mask=pd.read_csv(analysis_folder+savefilename[:-4]+"_filtermask.tsv",sep="\t",index_col=list(range(len(ids))),header=[0,1])
        mask_description=pd.read_csv(analysis_folder+savefilename[:-4]+"_filtermask_description.tsv",sep="\t").set_index("filter")
        print("Existing profiles found and loaded.")
    except (FileNotFoundError, KeyError):
        print("Loading raw DIA-NN output, filtering and creating profiles")
        if filters == dict():
            print("No 'filters' dictionary provided, analysis run without any filtering")
            filters = {"filter_proteotypic": False, "filter_by_qvalue": False, "filter_by_global_qvalue": False, "filter_by_gg_qvalue": False, "filter_by_pg_qvalue": False, "filter_by_p_qvalue": False, "filter_by_lib_qvalue": False, "filter_by_PEP": False, "filter_by_MScorr": False,
                       "filter_by_evidence": False, "filter_by_MScorr_evidence_combined": False, "qvalue_thresh": 0.01, "lib_qvalue_thresh": 0.01, "pep_thresh": 0.05, "MScorr_thresh": 0.9, "evidence_thresh": 1, "MScorr_evidence_combined_thresh": 0.9, "min_peptides": 2}
        if profile_annotation_file is None:
            profile_annotation_file = savefilename[:-4] + \
                "_filename_annotations.tsv"
        # load raw file
        file, path, filename = read_input(filename=raw_filename)
        #file = filtering(file=file, filter_proteotypic=filters["filter_proteotypic"], filter_by_qvalue=filters["filter_by_qvalue"], filter_by_global_qvalue=filters["filter_by_global_qvalue"], filter_by_gg_qvalue=filters["filter_by_gg_qvalue"], filter_by_pg_qvalue=filters["filter_by_pg_qvalue"], filter_by_p_qvalue=filters["filter_by_p_qvalue"], filter_by_lib_qvalue=filters["filter_by_lib_qvalue"], filter_by_PEP=filters["filter_by_PEP"],filter_by_MScorr=filters["filter_by_MScorr"], filter_by_evidence=filters["filter_by_evidence"], filter_by_MScorr_evidence_combined=filters["filter_by_MScorr_evidence_combined"], qvalue_thresh=filters["qvalue_thresh"], lib_qvalue_thresh=filters["lib_qvalue_thresh"], pep_thresh=filters["pep_thresh"], MScorr_thresh=filters["MScorr_thresh"], evidence_thresh=filters["evidence_thresh"], MScorr_evidence_combined_thresh=filters["MScorr_evidence_combined_thresh"])
        module_helper.create_analysis_folder(analysis_folder)
        # filering on peptide level
        # file=filtering_peptide_specific(file=file) #UNDER CONSTRUCTION
        profile_annotation_df = module_helper.create_filename_annotations_file(path=path, file=file, profile_annotation_file=profile_annotation_file,
                                                                 profile_order=profile_order, raw_data_folder=raw_data_folder, profile_annotation=profile_annotation, filename2run_pattern=filename2run_pattern, filename2samplegroup_pattern=filename2samplegroup_pattern, overwrite_filename_annotations=overwrite_filename_annotations,encoding=encoding, diann_sample_column=diann_sample_column)
        print("If wanted, annotate names of runs in "+analysis_folder+profile_annotation_file +
              " in column name on x-axis.\n You can reorder the runs as well.")
        input("Press enter when finished annotations (saved and closed file) or to skip...")
        profile_annotation_df = module_helper.create_filename_annotations_file(path=path, file=file, profile_annotation_file=profile_annotation_file,
                                                                 profile_order=profile_order, raw_data_folder=raw_data_folder, profile_annotation=profile_annotation, filename2run_pattern=filename2run_pattern, filename2samplegroup_pattern=filename2samplegroup_pattern, overwrite_filename_annotations=overwrite_filename_annotations,encoding=encoding,diann_sample_column=diann_sample_column)
        file,filter_mask,mask_description = filtering(file=file, filters=filters, ids=ids, profile_annotation_df=profile_annotation_df, values=values, filename_col=filename_col, diann_sample_column=diann_sample_column)
        profiles_file, profiles_file_numpeptides, filter_mask, mask_description = create_profiles(file=file, ids=ids, values=values, profile_annotation_df=profile_annotation_df,
                                                                   analysis_folder=analysis_folder, savefilename=savefilename[:-4], filter_mask=filter_mask,mask_description=mask_description,min_peptides=filters["min_peptides"], count_peptides=filters["count_peptides"], diann_sample_column=diann_sample_column)
    return([profiles_file, profiles_file_numpeptides, filter_mask, mask_description])



