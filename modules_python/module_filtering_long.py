# -*- coding: utf-8 -*-"
"""
Input: Data import from software e.g. DIA-NN
Output: Filtered wide-format data

         
         
functions:
    - read_input: reads DIA-NN output report.tsv or MSFragger output
    - select_filtering_columns: provides based on applied filters a pre-selection of column names to read from the DIA-NN output report.tsv to save RAM
    - filtering: apply filtering on DIA-NN output report.tsv
    (- filtering_peptide_specific: some optional additional filters based on peptide properties (so far only based on precursor charge state) => not used)
    - name_unknown_proteins: replaces NaN in id_mapping_df[protein_col], i.e. renames missing protein origins.
    - clean_mask_description_df: removes filters which were not applied from mask_description_df
    (- create_profiles: converts (filtered) DIA-NN output report.tsv to wide-format profiles with IDs as rows and samples as columns (as specified by metadata) => Can be optimised in runtime)
    - long_to_wide _format: Transforms long format data frame to wide format


variable names + types
    DataFrame: *_df
    array: *_arr
    dictionary: *_dict
    list: *_list
    set: *_set
    int/float/bool/str: *
    filenames: *_fp

#TODO
include new possible filters by DIA-NN >= 2.

"""
# %%


import sys
import pandas as pd
import numpy as np
import time
import os
import re
#from modules_python import module_helper as helper
import module_helper as helper



##############################################################################
# File handling
##############################################################################


def read_input(filename,file_columns_list=None):
    """
    Reads DIA-NN output report.tsv or MSFragger output, if specified only provided columns.
    author: BZ
    
    Parameters
    ----------
    filename : str, optional
        string specifying DIA-NN output path and filename (report.tsv). The default is None, opening a GUI to select the file.
    file_columns_list : list, optional
        list specifying which columns to read from DIA-NN output. The default is None, reading all columns.

    Returns
    -------
    pandas.DataFrame

    """
        
    #filesize = os.path.getsize(filename)
    #ram=psutil.virtual_memory()
    #if filesize<ram.available:
    #read parquet file
    if "parquet" in filename.split("."):
         file = pd.read_parquet(filename, columns=file_columns_list)
    else: # read tsv file
         file = pd.read_csv(filename, usecols=file_columns_list, sep="\t")# ,dtype={"Q.Value":float,"GG.Q.Value":float,"PG.Q.Value":float,"Protein.Q.Value":float,"Lib.Q.Value":float})
    
    return(file)




def select_filtering_columns(filters_dict,ids_list=["Precursor.Id", "Stripped.Sequence", "Protein.Group", "Genes"],values_col="Precursor.Quantity",input_filename_col="File.Name"):
    """
    with given filters dictionary, ids and values columns this function creates a list with column names of interest for extraction from DIA-NN output report.tsv
    author: BZ
    
    Parameters
    ----------
    filters_dict : dict
        dictionary containing the used filters and filtering parameters.
    ids_list : list
        list of ids extracted from DIA-NN report.tsv. The default is ["Precursor.Id", "Stripped.Sequence", "Protein.Group", "Genes"].
    values_col : str
        values extracted from DIA-NN report.tsv, e.g. Precursor.Quantity.
    input_filename_col : str, optional
        name of the sample names in the DIA-NN output, DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    
    Returns
    -------
    file_columns_list : list
        list of e.g. DIANN output report.tsv column names to read.

    """
    #Only read necessary columns to save RAM (these computers are weak)
    filtering_columns_list=list()
    if filters_dict["filter_proteotypic"]:
        filtering_columns_list.append("Proteotypic")
    if filters_dict["filter_by_qvalue"]:
        filtering_columns_list.append("Q.Value")
    if filters_dict["filter_by_global_qvalue"]:
        filtering_columns_list.append("Global.Q.Value")
    if filters_dict["filter_by_gg_qvalue"]:
        filtering_columns_list.append("GG.Q.Value")
    if filters_dict["filter_by_pg_qvalue"]:
        filtering_columns_list.append("PG.Q.Value")
    if filters_dict["filter_by_global_pg_qvalue"]:
        filtering_columns_list.append("Global.PG.Q.Value")
    if filters_dict["filter_by_p_qvalue"]:
        filtering_columns_list.append("Protein.Q.Value")
    if filters_dict["filter_by_lib_qvalue"]:
        filtering_columns_list.append("Lib.Q.Value")
    if filters_dict["filter_by_lib_pg_qvalue"]:
        filtering_columns_list.append("Lib.PG.Q.Value")
    if filters_dict["filter_by_PEP"]:
        filtering_columns_list.append("PEP")
    if filters_dict["filter_by_MScorr"] or filters_dict["filter_by_MScorr_evidence_combined"]:
        filtering_columns_list.append("Ms1.Profile.Corr")
    if filters_dict["filter_by_evidence"] or filters_dict["filter_by_MScorr_evidence_combined"]:
        filtering_columns_list.append("Mass.Evidence")
    file_columns_list = [input_filename_col] + ids_list + [values_col] + filtering_columns_list
    return(file_columns_list)




##############################################################################
# FILTERING
##############################################################################

# filter DIA-NN output (except on number of peptides), by default all filters set to False (no filtering)
def filtering(input_df, filters_dict, input_id_col="Precursor.Id", metadata_df=None, mask_df=None, mask_description_df=None, keep_track=True, values_col="Precursor.Quantity", filename_col="filename", diann_sample_column="File.Name"):
    """
    Filters rows in DIA-NN output report.tsv on filters as specified in filters_dict by True.
    Number of filtered precursors in every step are saved to logfile.
    
    author: BZ
    
    Parameters
    ----------
    input_df : pd.DataFrame
        E.g. DIA-NN output report.tsv, long format.
    filters_dict : dict
        dictionary containing the used filters and filtering parameters.
    input_id_col : str
        name of column with IDs extracted from DIA-NN report.tsv. The default is "Precursor.Id".
    metadata_df : pd.DataFrame
        Metadata mapping File.Name to shorter run names. The default is None. 
    mask_df : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, calling creating_filtering_matrix() to create this variable.
    mask_description_df: pd.DataFrame.
        mapping numbers in mask_df to filters. The default is None, calling creating_filtering_matrix() to create this variable.
    keep_track : bool, optional
        Flag to activate/deactivate tracking of filtered values in mask_df. The default is True.
    values_col : str
        values extracted from DIA-NN report.tsv. The default is "Precursor.Quantity".
    filename_col : str
        name of the column in metadata_df specifying filenames (including path) as in DIA-NN report.tsv. The default is "filename".
    diann_sample_column : str, optional
        name of the sample names in the DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    Returns
    -------
    input_df : pd.DataFrame
        filtered DIA-NN output report.tsv.
    mask_df : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.
    mask_description_df : pd.DataFrame.
        mapping numbers in mask_df to filters.
    log : str
        string containing messages for report log.

    """
    log = ""
    len_input_df_before = len(input_df)
    #Exclude runs from input_df not present in metadata
    if metadata_df is not None:
        intersect_files = input_df[diann_sample_column].isin(metadata_df[filename_col])
        if intersect_files.sum() == 0:
            intersect_files = input_df[diann_sample_column].isin(metadata_df.index)
        if intersect_files.sum() == 0:
            log += "\t\t\t\t###WARNING###: No overlap between samples in metadata and DIA-NN report. please check wheter filename_col is defined correctly."
            #sys.exit("###WARNING###: No overlap between samples in metadata and DIA-NN report. please check wheter filename_col is defined correctly.")
            return([None,None,None,log])
        input_df = input_df[intersect_files]
        if len_input_df_before > len(input_df):
            log += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t###WARNING###: Not all files in input_df are found in metadata. Analysis performed only on files present in metadata.\n"
        if keep_track:
            if mask_df is None or mask_description_df is None:
                mask_df,mask_description_df,log_create = helper.create_filtering_matrix(input_df,metadata_df,min_peptides=filters_dict["min_peptides"],min_precursors=filters_dict["min_precursors"],values=values_col,ids=input_id_col,filename_col=filename_col,diann_sample_column=diann_sample_column)
                log += log_create
    # FILTERING
    filtering_tracker_df = pd.DataFrame(index={key:value for key,value in filters_dict.items() if value==True}.keys(),columns=["number of excluded precursors","number of remaining precursors"])
    len_input_df_before = len(input_df)
    len_input_df_after = len(input_df)
    log += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tBegin Filtering\n"
    # FILTER BY PROTEOTYPIC
    # exclude non-proteotypic
    if filters_dict["filter_proteotypic"]:
        if "Proteotypic" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["Proteotypic"]!=1][input_id_col].tolist()
                cols = input_df[input_df["Proteotypic"]!=1][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]]=mask_description_df.loc["Proteotypic","value"]
            input_df = input_df[input_df["Proteotypic"] == 1]
            filtering_tracker_df.loc["filter_proteotypic","number of excluded precursors"] = len_input_df_after-len(input_df)
            filtering_tracker_df.loc["filter_proteotypic","number of remaining precursors"] = len(input_df)
            len_input_df_after=len(input_df)
        else:
            log += "\t\t\t\tColumn 'Proteotypic' not in input_df, filter skipped\n"
    # FILTER BY Q-VALUE
    # exclude by q-value threshold for gene group AND precursor
    if filters_dict["filter_by_qvalue"]:
        if "Q.Value" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["Q.Value"] > filters_dict["qvalue_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["Q.Value"] > filters_dict["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["Q.Value","value"]
            input_df = input_df[input_df["Q.Value"] <= filters_dict["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_q_value","number of excluded precursors"] = len_input_df_after-len(input_df)
            filtering_tracker_df.loc["filter_by_q_value","number of remaining precursors"] = len(input_df)
            len_input_df_after=len(input_df)
        else:
            log += "\t\t\t\tColumn 'Q.Value' not in input_df, filter skipped\n"
    if filters_dict["filter_by_global_qvalue"]:
        if "Global.Q.Value" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["Global.Q.Value"] > filters_dict["qvalue_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["Global.Q.Value"]>filters_dict["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["Global.Q.Value","value"]
            input_df = input_df[input_df["Global.Q.Value"] <= filters_dict["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_global_qvalue","number of excluded precursors"] = len_input_df_after-len(input_df)
            filtering_tracker_df.loc["filter_by_global_qvalue","number of remaining precursors"] = len(input_df)
            len_input_df_after = len(input_df)
        else:
            log += "\t\t\t\tColumn 'Global.Q.Value' not in input_df, filter skipped\n"
    if filters_dict["filter_by_gg_qvalue"]:
        if "GG.Q.Value" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["GG.Q.Value"]>filters_dict["qvalue_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["GG.Q.Value"]>filters_dict["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["GG.Q.Value","value"]
            input_df = input_df[input_df["GG.Q.Value"] <= filters_dict["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_gg_qvalue","number of excluded precursors"]=len_input_df_after-len(input_df)
            filtering_tracker_df.loc["filter_by_gg_qvalue","number of remaining precursors"]=len(input_df)
            len_input_df_after=len(input_df)
        else:
            log += "\t\t\t\tColumn 'GG.Q.Value' not in input_df, filter skipped\n"
    if filters_dict["filter_by_pg_qvalue"]:
        if "PG.Q.Value" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["PG.Q.Value"] > filters_dict["qvalue_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["PG.Q.Value"] > filters_dict["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["PG.Q.Value","value"]
            input_df = input_df[input_df["PG.Q.Value"] <= filters_dict["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_pg_qvalue","number of excluded precursors"] = len_input_df_after - len(input_df)
            filtering_tracker_df.loc["filter_by_pg_qvalue","number of remaining precursors"] = len(input_df)
            len_input_df_after=len(input_df)
        else:
            log += "\t\t\t\tColumn 'PG.Q.Value' not in input_df, filter skipped\n"
    if filters_dict["filter_by_p_qvalue"]:
        if "Protein.Q.Value" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["Protein.Q.Value"] > filters_dict["qvalue_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["Protein.Q.Value"] > filters_dict["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["Protein.Q.Value","value"]
            input_df = input_df[input_df["Protein.Q.Value"] <= filters_dict["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_p_qvalue","number of excluded precursors"] = len_input_df_after - len(input_df)
            filtering_tracker_df.loc["filter_by_p_qvalue","number of remaining precursors"] = len(input_df)
            len_input_df_after = len(input_df)
        else:
            log += "\t\t\t\tColumn 'Protein.Q.Value' not in input_df, filter skipped\n"
    if filters_dict["filter_by_global_pg_qvalue"]:
        if "Global.PG.Q.Value" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["Global.PG.Q.Value"] > filters_dict["qvalue_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["Global.PG.Q.Value"] > filters_dict["qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]]=mask_description_df.loc["Global.PG.Q.Value","value"]
            input_df = input_df[input_df["Global.PG.Q.Value"] <= filters_dict["qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_global_pg_qvalue","number of excluded precursors"]=len_input_df_after - len(input_df)
            filtering_tracker_df.loc["filter_by_global_pg_qvalue","number of remaining precursors"]=len(input_df)
            len_input_df_after=len(input_df)
        else:
            log += "\t\t\t\tColumn 'Global.PG.Q.Value' not in input_df, filter skipped\n"
    if filters_dict["filter_by_lib_qvalue"]:
        if "Lib.Q.Value" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["Lib.Q.Value"] > filters_dict["lib_qvalue_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["Lib.Q.Value"] > filters_dict["lib_qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["Lib.Q.Value","value"]
            input_df = input_df[input_df["Lib.Q.Value"] <= filters_dict["lib_qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_lib_qvalue","number of excluded precursors"] = len_input_df_after - len(input_df)
            filtering_tracker_df.loc["filter_by_lib_qvalue","number of remaining precursors"] = len(input_df)
            len_input_df_after = len(input_df)
        else:
            log += "\t\t\tColumn 'Lib.Q.Value' not in input_df, filter skipped\n"
    if filters_dict["filter_by_lib_pg_qvalue"]:
        if "Lib.PG.Q.Value" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["Lib.PG.Q.Value"] > filters_dict["lib_qvalue_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["Lib.PG.Q.Value"] > filters_dict["lib_qvalue_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]]=mask_description_df.loc["Lib.PG.Q.Value","value"]
            input_df = input_df[input_df["Lib.PG.Q.Value"] <= filters_dict["lib_qvalue_thresh"]]
            filtering_tracker_df.loc["filter_by_lib_pg_qvalue","number of excluded precursors"] = len_input_df_after - len(input_df)
            filtering_tracker_df.loc["filter_by_lib_pg_qvalue","number of remaining precursors"] = len(input_df)
            len_input_df_after = len(input_df)
        else:
            log += "\t\t\t\tColumn 'Lib.PG.Q.Value' not in input_df, filter skipped\n"

    # FILTER BY PEP
    # exclude rows with posterior error probability (PEP) for precursor identification exceeding pep_thresh
    if filters_dict["filter_by_PEP"]:  # NOT CHECKED
        if "PEP" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["PEP"] > filters_dict["pep_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["PEP"] > filters_dict["pep_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["PEP","value"]
            input_df = input_df[input_df["PEP"] <= filters_dict["pep_thresh"]]
            filtering_tracker_df.loc["filter_by_PEP","number of excluded precursors"] = len_input_df_after - len(input_df)
            filtering_tracker_df.loc["filter_by_PEP","number of remaining precursors"] = len(input_df)
            len_input_df_after = len(input_df)
        else:
            log += "\t\t\t\tColumn 'PEP' not in input_df, filter skipped\n"
    # FILTER ON MS-LEVEL
    # exclude rows with correlation between MS1 and MS2 chromatograms below MScorr_thresh
    if filters_dict["filter_by_MScorr"]:
        if "Ms1.Profile.Corr" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["Ms1.Profile.Corr"] < filters_dict["MScorr_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["Ms1.Profile.Corr"] < filters_dict["MScorr_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["Ms1.Profile.Corr","value"]
            input_df = input_df[input_df["Ms1.Profile.Corr"] >= filters_dict["MScorr_thresh"]]
            filtering_tracker_df.loc["filter_by_MScorr","number of excluded precursors"] = len_input_df_after - len(input_df)
            filtering_tracker_df.loc["filter_by_MScorr","number of remaining precursors"] = len(input_df)
            len_input_df_after = len(input_df)
        else:
            log += "\t\t\t\tMs1.Profile.Corr not in input_df, filter skipped\n"
    # exclude rows with MS2 evidence below evidence_thresh
    if filters_dict["filter_by_evidence"]:
        if "Mass.Evidence" in input_df.columns:
            if keep_track:
                rows = input_df[input_df["Mass.Evidence"] < filters_dict["evidence_thresh"]][input_id_col].tolist()
                cols = input_df[input_df["Mass.Evidence"] < filters_dict["evidence_thresh"]][diann_sample_column].tolist()
                for r in range(len(rows)):
                    mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["Mass.Evidence","value"]
            input_df = input_df[input_df["Mass.Evidence"] >= filters_dict["evidence_thresh"]]
            filtering_tracker_df.loc["filter_by_evidence","number of excluded precursors"] = len_input_df_after-len(input_df)
            filtering_tracker_df.loc["filter_by_evidence","number of remaining precursors"] = len(input_df)
            len_input_df_after = len(input_df)
        else:
            log += "\t\t\t\tMass.Evidence not in input_df, filter skipped\n"
    if filters_dict["filter_by_MScorr_evidence_combined"]:
        if "Ms1.Profile.Corr" in input_df.columns and "Mass.Evidence" in input_df.columns:
            if keep_track:
                log += "###WARNING###: tracking filtering for filter_by_MScorr_evidene_combined not implemented yet"#IMPLEMENT IF NECCESSARY
            input_df = input_df[input_df["Ms1.Profile.Corr"] + 0.5 *
                        input_df["Mass.Evidence"] >= filters_dict["MScorr_evidence_combined_thresh"]]
            filtering_tracker_df.loc["filter_by_MScorr_evidence_combined","number of excluded precursors"] = len_input_df_after - len(input_df)
            filtering_tracker_df.loc["filter_by_MScorr_evidence_combined","number of remaining precursors"] = len(input_df)
            len_input_df_after = len(input_df)
        else:
            log += "\t\t\t\tMs1.Profile.Corr and/or Mass.Evidence not in input_df, filter skipped\n"
    #filter out precursors with zero value
    if filters_dict["filter_by_zero_values"]:
        if keep_track:
            rows = input_df[input_df[values_col] == 0][input_id_col].tolist()
            cols = input_df[input_df[values_col] == 0][diann_sample_column].tolist()
            for r in range(len(rows)):
                mask_df.loc[rows[r],cols[r]] = mask_description_df.loc["Zero "+values_col,"value"]
        input_df=input_df[input_df[values_col] > 0]
        filtering_tracker_df.loc["filter_by_zero_values","number of excluded precursors"] = len_input_df_after - len(input_df)
        filtering_tracker_df.loc["filter_by_zero_values","number of remaining precursors"] = len(input_df)
        len_input_df_after = len(input_df)
    len_input_df_after = len(input_df)
    filtering_tracker_df.dropna(inplace=True)
    log += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tFiltering completed. "+str(len_input_df_before-len_input_df_after)+"("+str(np.round(100-len_input_df_after/len_input_df_before*100,2))+"%) Peptides were excluded during filtering of input_df\n\n"
    log += filtering_tracker_df.to_string().replace("\n","\n\t") + "\n"
    if input_df.index.tolist() == []:
        log += "\t\t\t\t###WARNING### Filtered input_df is empty. Please check if columns on which filters are based are not empty (sometimes depending on library used e.g. Genes and therefore Q-Values on gene-level are missing), else adjust filtering parameters.\n"
    return([input_df,mask_df,mask_description_df,log])



#Construction site, improve if needed
'''
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
'''

#Not Tested #TODO
def name_unknown_proteins(id_mapping_df,input_fp, ids_list=["Precursor.Id","Stripped.Sequence","Protein.Group","Genes"], protein_col="Protein.Group", replace_id_col="Protein.Ids"):
    """
    Replaces NaN in protein_col of id_mapping_df with values from the column replace_id_col in the dataframe from which id_mapping was initially generated (file at input_fp).
    Remaining NaN are replaced with "UnknownProteinX" where X is a unique integer, thus loosing the information whether Precursors with missing protein_col belonged to the same protein.
    author: BZ
    
    Parameters
    ----------
    id_mapping_df : pd.DataFrame
        DataFrame collecting IDs such as Gene Names or Protein Groups for Precursors
    input_fp : str
        path and filename of a file mapping Precursors in id_mapping_df to other IDs, e.g. the DIA-NN output report.tsv.
    ids_list : list, optional
        List of columns in id_mapping_df used to uniquely identify a precursor. The default is ["Precursor.Id","Stripped.Sequence","Protein.Group","Genes"].
    protein_col : str, optional
        Name of the column in id_mapping_df and the data stored at input_fp with missing IDs which are to be replaced. The default is "Protein.Group".
    replace_id_col : str, optional
        Name of the column in the data stored at input_fp to be used for replacement of missing ids in id_mapping_df[protein_col]. The default is "Protein.Ids".

    Returns
    -------
    id_mapping_df : pd.DataFrame
        DataFrame collecting IDs such as Gene Names or Protein Groups for Precursors.

    """
    #define columns needed from input
    ids_list = list(set([protein_col] + ids_list))
    file_columns = list(set([replace_id_col] + ids_list))
    #read input to find possible id translation
    trans_df = read_input(input_fp, file_columns=file_columns).set_index(ids_list)
    #replace NaN in protein_col with Protein.Ids from input
    protein_groups = [trans_df.loc[x,replace_id_col] if pd.isna(x[protein_col]) else x for x in id_mapping_df[ids_list]]
    #if there are still NaN, replace them with UnknownProtein + some integer, 
    i=1
    for x in range(len(protein_groups)):
        if pd.isna(protein_groups[x]):
            protein_groups[x] = "UnknownProtein"+str(i)
            i+=1
    #assign new values to initial mapping
    id_mapping_df[protein_col] = protein_groups
    return(id_mapping_df)

def create_peptide_ids(precursors):
    """
    strips precursors from charge state to create peptide ids.
    author: BZ
    
    Parameters
    ----------
    precursors : list
        list of precursor ids.

    Returns
    -------
    list
        peptide ids matching precursors.

    """
    return [re.sub(r"\d+$","",x) for x in precursors]


def clean_mask_description_df(mask_description_df,filters_dict):
    """
    removes filters which were not applied from mask_description_df.
    author: BZ

    Parameters
    ----------
    mask_description_df: pd.DataFrame.
        mapping numbers in filtering mask to filters.
    filters_dict : dict
        dictionary containing the used filters and filtering parameters.

    Returns
    -------
    mask_description_df: pd.DataFrame.
        mapping numbers in mask to filters.
        
    """
    if not filters_dict["filter_proteotypic"]:
        mask_description_df.drop("Proteotypic",axis=0,inplace=True)
    if not filters_dict["filter_by_qvalue"]:
        mask_description_df.drop("Q.Value",axis=0,inplace=True)
    if not filters_dict["filter_by_global_qvalue"]:
        mask_description_df.drop("Global.Q.Value",axis=0,inplace=True)
    if not filters_dict["filter_by_gg_qvalue"]:
        mask_description_df.drop("GG.Q.Value",axis=0,inplace=True)
    if not filters_dict["filter_by_pg_qvalue"]:
        mask_description_df.drop("PG.Q.Value",axis=0,inplace=True)
    if not filters_dict["filter_by_global_pg_qvalue"]:
        mask_description_df.drop("Global.PG.Q.Value",axis=0,inplace=True)
    if not filters_dict["filter_by_p_qvalue"]:
        mask_description_df.drop("Protein.Q.Value",axis=0,inplace=True)
    if not filters_dict["filter_by_lib_qvalue"]:
        mask_description_df.drop("Lib.Q.Value",axis=0,inplace=True)
    if not filters_dict["filter_by_lib_pg_qvalue"]:
        mask_description_df.drop("Lib.PG.Q.Value",axis=0,inplace=True)
    if not filters_dict["filter_by_PEP"]:
        mask_description_df.drop("PEP",axis=0,inplace=True)
    if not filters_dict["filter_by_MScorr"]:
        mask_description_df.drop("Ms1.Profile.Corr",axis=0,inplace=True)
    if not filters_dict["filter_by_evidence"]:
        mask_description_df.drop("Mass.Evidence",axis=0,inplace=True)    
    if not filters_dict["filter_by_MScorr_evidence_combined"]:
        mask_description_df.drop("Ms1.Profile.Corr + Mass.Evidence",axis=0,inplace=True)
        
    return mask_description_df


#DEPRECATED
'''
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
        filter_mask,mask_description=helper.create_filtering_matrix(file,profile_annotation_df,min_peptides=min_peptides,values=values,ids=ids,diann_sample_column=diann_sample_column)
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
        print("###WARNING###: Multiple times occuring IDs found (ok on protein level)")
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
            
            #ID=index[1]
            #xvalue=index[0]
            #if np.round(file.loc[(xvalue, ID), values][0],0)!=np.round(file.loc[(xvalue, ID), values].mean(),0): 
            #    print("oh-oh, different values over peptides encountered")
            
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
'''



def long_to_wide_format(long_df,metadata_df,mask_df,mask_description_df,id_col="Precursor.Id",value_col="Precursor.Quantity",input_filename_col="File.Name",metadata_filename_col="File.Name",metadata_sample_col="Sample.Name"):
    """
    Transforms long format data frame to wide format. 
    author: BZ

    Parameters
    ----------
    long_df : pd.DataFrame
        (filtered) DIA-NN output report.tsv.
    metadata_df : pd.DataFrame
        etadata mapping diann_sample_column (File.Name/Run) to shorter run names.
    mask_df : pd.DataFrame
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.
    mask_description_df : pd.DataFrame
        mapping numbers in mask to filters
    id_col : str, optional
        column name in long_df with precursor ids. The default is "Precursor.Id".
    value_col : str, optional
        column name in long_df with quantities. The default is "Precursor.Quantity".
    input_filename_col : str, optional
        column name in long_df with file/sample names. The default is "File.Name".
    metadata_filename_col : str, optional
        column name with complete filenames in metadata, which match file names in long_df[input_filename_col]. The default is "File.Name".
    metadata_sample_col : str, optional
        column name with (shorter) sample names in metadata. The default is "Sample.Name".

    Returns
    -------
    wide_df : pd.DataFrame
        wide format dataframe, rows contain IDs as described by id_col, columns contain samples as described by metdata_sample_col and values are values from the column value_col in long_df.
    mask_df : pd.DataFrame
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered and which were measured.
    """
    #Sort index
    long_df.sort_index(inplace=True)
    #transform to wide format
    wide_df = pd.pivot_table(long_df.reset_index(),index=id_col,columns=input_filename_col,values=value_col)
    #rename columns from file names to sample names
    wide_df.columns = metadata_df.reset_index().set_index(metadata_filename_col,drop=False).loc[wide_df.columns,metadata_sample_col]
    wide_df.sort_index(inplace=True,axis=1)
    #assign measured values in mask_df
    for col in wide_df.columns:
        mask_df.loc[wide_df[wide_df[col]>0].index,col] = mask_description_df.loc["measured and kept","value"]
    return([wide_df,mask_df])
    
    
# %%


if __name__ == "__main__":
    
    #temporary for debugging while load_parameters is not ready, careful with eval if variables are already defined, quick and dirty
    '''
    out_dir_fp = "../01_Filtering/"
    parameters_dict_raw = pd.read_csv("../01_Filtering/parameters.txt",sep="\t",index_col=0,header=None)[1].to_dict()
    parameters_dict = dict()
    for x in parameters_dict_raw.items():
        try: #if type != str
            parameters_dict[x[0]] = eval(x[1])
        except (SyntaxError, NameError): #if type==str
            parameters_dict[x[0]] = x[1]
    
    data_df_fp = parameters_dict["input_fp"]
    metadata_df,log = helper.load_metadata(metadata_fp=parameters_dict["metadata_fp"],
                                              metadata_sheet=parameters_dict["metadata_sheet"], 
                                              filename_col=parameters_dict["metadata_filename_column"], 
                                              encoding=parameters_dict["metadata_encoding"])
 
    '''
    
    data_df_fp, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, mask_description_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = "Filtering long", 
                                                                                                                                                                              load_metadata=True, 
                                                                                                                                                                              load_id_mapping=False, 
                                                                                                                                                                              load_mask=False, 
                                                                                                                                                                              return_data_fp_only=True)
    
    # Check if loading was successful, else pass error to log and exit   
    helper.logger(log_file,log)
    if metadata_df is None:
        sys.exit("###ERROR###: Metadata loading failed.")
    if "###WARNING###" in log:
        helper.logger(log_file,"###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])
        sys.exit("###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])

    # Log parameters
    log = time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Filtering Long Format\n"
    log += "\nParameters\n------------\n"
    log += "\t" + pd.DataFrame(parameters_dict.items(),columns=["Parameter","Value"]).to_string(index=False).replace("\n","\n\t") 
    log += "\n------------\n"  
    helper.logger(log_file, module_output=log)
         
    # Extract variables from parameters, which need adjustment
    ids_list = parameters_dict["input_ids_columns"].split(";")
    input_id_col = ids_list[0]
    
   
    
    #DIA-NN
    if parameters_dict["input_origin"]=="diann":
        # Select columns which will be read from DIA-NN output (reading only used columns saves RAM and speeds up analysis)
        file_columns_list = select_filtering_columns(filters_dict=parameters_dict, 
                                                     ids_list=ids_list, 
                                                     values_col=parameters_dict["input_values_column"],
                                                     input_filename_col=parameters_dict["input_filename_column"])
    
        
        
        try:
            #Loading input file
            log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tLoading input file " + data_df_fp + "\n"
            helper.logger(log_file,log)
            input_df = read_input(filename=data_df_fp, file_columns_list=file_columns_list)
        except KeyError:
            #read complete input
            input_df = read_input(filename=data_df_fp, file_columns_list=None)
            #identify missing columns which lead to KeyError; DIA-NN specific!
            missing_cols_set = set(file_columns_list) - set(input_df.columns)
            columns_filter_dict = {"Proteotypic":"filter_proteotypic",
                                   "Q.Value":"filter_by_qvalue",
                                   "Global.Q.Value":"filter_by_global_qvalue",
                                   "GG.Q.Value":"filter_by_gg_qvalue","PG.Q.Value":"filter_by_pg_qvalue",
                                   "Global.PG.Q.Value":"filter_by_global_pg_qvalue",
                                   "Protein.Q.Value":"filter_by_p_qvalue",
                                   "Lib.Q.Value":"filter_by_lib_qvalue",
                                   "Lib.PG.Q.Value":"filter_by_lib_pg_qvalue",
                                   "PEP":"filter_by_PEP",
                                   "Ms1.Profile.Corr":"filter_by_MScorr",
                                   "Mass.Evidence":"filter_by_evidence"}
            if missing_cols_set - set(columns_filter_dict.keys()) == set():
                #switch off filters linked to missing columns
                for missing_col in missing_cols_set:
                    parameters_dict[columns_filter_dict[missing_col]] = False
                    log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "###WARNING###"+"Filter " + columns_filter_dict[missing_col] + " set to False, because the column(s) "+", ".join(list(missing_cols_set) + " is/are missing in the DIA-NN output file.")
            if len(missing_cols_set & set(["Mass.Evidence","Ms1.Profile.Corr"])) > 0:
                parameters_dict["filter_by_MScorr_evidence_combined"]=False
            else:
                #if missing columns are linked to ids, values, or sample columns raise a KeyError and exit
                log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\t###ERROR###"+", ".join(missing_cols_set) + " not in the selected input file.\n"
                helper.logger(log_file,log)
                #raise KeyError(", ".join(missing_cols_set) + " not in the selected input file.")
                sys.exit("###ERROR###" + ", ".join(missing_cols_set) + " not in the selected input file.")
        log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tFile loaded successfully.\n"
        helper.logger(log_file,log)
        
        
        # Extract ID mapping
        id_mapping_df = input_df[ids_list + ["Proteotypic"]].drop_duplicates()

        # Rename missing values in protein_col
        if pd.isna(id_mapping_df[parameters_dict["protein_col"]]).sum() > 0:
            id_mapping_df = name_unknown_proteins(id_mapping_df=id_mapping_df,
                                                  input_fp=data_df_fp, 
                                                  ids_list=ids_list, 
                                                  protein_col=parameters_dict["protein_col"], 
                                                  replace_id_col=parameters_dict["alternative_protein_col"])
        # Create a column "Peptide.Id" with peptide ids in id_mapping_df, which contains precursor ids without charge state
        if parameters_dict["peptide_col"] not in id_mapping_df.columns:
            id_mapping_df[parameters_dict["peptide_col"]] = create_peptide_ids(id_mapping_df[parameters_dict["precursor_col"]])
                 
        # Precursor filtering (DIA-NN specific)
        input_df, mask_df, mask_description_df, log = filtering(input_df=input_df, 
                                                                  input_id_col=input_id_col, 
                                                                  metadata_df=metadata_df, 
                                                                  filters_dict=parameters_dict, 
                                                                  values_col=parameters_dict["input_values_column"], 
                                                                  filename_col=parameters_dict["metadata_filename_column"], 
                                                                  diann_sample_column=parameters_dict["input_filename_column"])
    elif parameters_dict["input_origin"]=="msfragger":
        pass # adjust here for other programs
        
    # Check if filtering was successful, else pass error to log and exit    
    if not input_df.empty:
        log += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tLong format filtering successful.\n"
        helper.logger(log_file,log)
    else:
        helper.logger(log_file,log)
        sys.exit("###ERROR###: Long format Filtering failed.")
        
    # filter mask_description for applied filters only
    mask_description_df = clean_mask_description_df(mask_description_df,filters_dict=parameters_dict)

    # Rename filter_mask columns
    mask_df = helper.rename_filter_mask_columns(mask_df, 
                                                metadata_df, 
                                                filename_col=parameters_dict["metadata_filename_column"], 
                                                sample_col=parameters_dict["metadata_sample_column"])
    
    # Long to wide format
    wide_df, mask_df = long_to_wide_format(input_df,
                                  metadata_df,
                                  mask_df, 
                                  mask_description_df,
                                  id_col=input_id_col,
                                  value_col=parameters_dict["input_values_column"],
                                  input_filename_col=parameters_dict["input_filename_column"],
                                  metadata_filename_col=parameters_dict["metadata_filename_column"],
                                  metadata_sample_col=parameters_dict["metadata_sample_column"])

    # Log Transformation (Base 2)
    wide_df = helper.log_transform(wide_df.replace(0,np.nan), base=2)
    log += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tWide format data logtransformed to base 2\n"
    helper.logger(log_file,log)
    
    # Save data
    id_mapping_df.to_csv(os.path.join(out_dir_fp, "id_mapping.tsv"), index=False, sep="\t")
    #wide_df.to_csv(out_dir_fp + "data_" + str(counter) + ".tsv", sep="\t")
    #mask_df.to_csv(out_dir_fp + "data_" + str(counter) + "_mask.tsv", sep="\t")
    #mask_description_df.to_csv(out_dir_fp + "data_" + str(counter) + "_mask_description.tsv", sep="\t")
    
    module_output = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tWide format matrix saved successfully to " + os.path.join(out_dir_fp, "data_" + str(counter) + ".tsv")+"\n\t\t\tmask saved successfully to " + os.path.join(out_dir_fp,"data_mask_" + str(counter) + ".tsv")+"\n\t\t\tand table for ID mapping saved to " + os.path.join(out_dir_fp, "id_mapping.tsv")+"\n"
    helper.close_module(data_processed_df=wide_df, masked_df=mask_df, mask_description_df=mask_description_df, log_file=log_file, module_output=module_output, out_dir_fp=out_dir_fp, outname=outname, mask=True)
  








