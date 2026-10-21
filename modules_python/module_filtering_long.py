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
    - long_to_wide_format: Transforms long format data frame to wide format


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
- Decoy == 0
- Peptidoform.Q.Value
- Lib.Peptidoform.Q.Value
- Global.Peptidoform.Q.Value
- Channel.Q.Value 0.01 - 0.5
- Quantity.Quality >= 0.5
- Empirical.Quality
- PG.MaxLFQ.Quality >) 0.7
- Normalisation.Factor
- Normalisation.Noise
- Genes.MaxLFQ.Quality
- Gene.MaxLFQ.Unique.Quality
- PTM.Site.Confidence
- Lib.PTM.Site.Confidence
- Translated.Q.Value
- PG.PEP
"""
# %%


import sys
import pandas as pd
import numpy as np
import time
import os
import re
import module_helper as helper
from scipy.sparse import coo_matrix
import operator

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
         file = pd.read_parquet(filename, columns=file_columns_list, engine="pyarrow")
    else: # read tsv file
         file = pd.read_csv(filename, usecols=file_columns_list, sep="\t")# ,dtype={"Q.Value":float,"GG.Q.Value":float,"PG.Q.Value":float,"Protein.Q.Value":float,"Lib.Q.Value":float})
    
    return(file)




def select_filtering_columns(filters_dict,values_col="Precursor.Quantity",input_filename_col="File.Name"):
    """
    with given filters dictionary, ids and values columns this function creates lists for filtering from DIA-NN output report.tsv
    author: BZ, OL
    
    Parameters
    ----------
    filters_dict : dict
        dictionary containing the used filters and filtering parameters.
    values_col : str
        values extracted from DIA-NN report.tsv, e.g. Precursor.Quantity.
    input_filename_col : str, optional
        name of the sample names in the DIA-NN output, DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    
    Returns
    -------
    filtering_operation_dict : dict
        dictionary containg all filters.
    """
    filtering_operation_dict = {}

    if filters_dict["filter_proteotypic"]:
        filtering_operation_dict.update({"filter_proteotypic": {
            "column": "Proteotypic",
            "operator": "==",
            "value": 1,
            "dtype": int}
        })
    
    if filters_dict["filter_by_qvalue"]:
        filtering_operation_dict.update({"filter_by_qvalue": {
            "column": "Q.Value",
            "operator": "<=",
            "value": filters_dict["qvalue_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_global_qvalue"]:
        filtering_operation_dict.update({"filter_by_global_qvalue": {
            "column": "Global.Q.Value",
            "operator": "<=",
            "value": filters_dict["qvalue_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_gg_qvalue"]:
        filtering_operation_dict.update({"filter_by_gg_qvalue": {
            "column": "GG.Q.Value",
            "operator": "<=",
            "value": filters_dict["qvalue_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_pg_qvalue"]:
        filtering_operation_dict.update({"filter_by_pg_qvalue": {
            "column": "PG.Q.Value",
            "operator": "<=",
            "value": filters_dict["qvalue_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_global_pg_qvalue"]:
        filtering_operation_dict.update({"filter_by_global_pg_qvalue": {
            "column": "Global.PG.Q.Value",
            "operator": "<=",
            "value": filters_dict["qvalue_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_p_qvalue"]:
        filtering_operation_dict.update({"filter_by_p_qvalue": {
            "column": "Protein.Q.Value",
            "operator": "<=",
            "value": filters_dict["qvalue_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_lib_qvalue"]:
        filtering_operation_dict.update({"filter_by_lib_qvalue": {
            "column": "Lib.Q.Value",
            "operator": "<=",
            "value": filters_dict["lib_qvalue_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_lib_pg_qvalue"]:
        filtering_operation_dict.update({"filter_by_lib_pg_qvalue": {
            "column": "Lib.PG.Q.Value",
            "operator": "<=",
            "value": filters_dict["lib_qvalue_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_PEP"]:
        filtering_operation_dict.update({"filter_by_PEP": {
            "column": "PEP",
            "operator": "<=",
            "value": filters_dict["pep_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_MScorr"] or filters_dict["filter_by_MScorr_evidence_combined"]:
        filtering_operation_dict.update({"filter_by_MScorr": {
            "column": "Ms1.Profile.Corr",
            "operator": ">=",
            "value": filters_dict["MScorr_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_evidence"] or filters_dict["filter_by_MScorr_evidence_combined"]:
        filtering_operation_dict.update({"filter_by_evidence": {
            "column": "Mass.Evidence",
            "operator": ">=",
            "value": filters_dict["evidence_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_MScorr_evidence_combined"]:
        filtering_operation_dict.update({"filter_by_MScorr_evidence_combined": {
            "column": "__dummy__mass__ms1__",
            "operator": ">=",
            "value": filters_dict["MScorr_evidence_combined_thresh"],
            "dtype": float}
        })
    
    if filters_dict["filter_by_zero_values"]:
        filtering_operation_dict.update({"filter_by_zero_values": {
            "column": values_col,
            "operator": ">",
            "value": 0,
            "dtype": float}
        })
        
        return(filtering_operation_dict)



            



##############################################################################
# FILTERING
##############################################################################

# DEPRECATED
'''
def exclude_runs_not_in_metadata(input_df, metadata_df=None, input_filename_column="File.Name", metadata_filename_col="filename"):
    """
    Trims input_df (input_df[diann_sample_column]) to same samples present in metadata_df[filename_col].
    
    author: BZ
    
    Parameters
    ----------
    input_df : pd.DataFrame
        E.g. DIA-NN output report.tsv, long format.
    metadata_df : pd.DataFrame
        Metadata mapping File.Name to shorter run names. The default is None. 
    input_filename_column : str, optional
        name of the sample names in the DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    metadata_filename_col : str
        name of the column in metadata_df specifying filenames (including path) as in DIA-NN report.tsv. The default is "filename".
    
    Returns
    -------
    input_df : pd.DataFrame
        filtered DIA-NN output report.tsv.
    log : str
        string containing messages for report log.
        
    """

    log = ""
    len_input_df_before = len(input_df)
    #Exclude runs from input_df not present in metadata
    if metadata_df is not None:
        intersect_files = input_df[input_filename_column].isin(metadata_df[metadata_filename_col])
        if intersect_files.sum() == 0:
            intersect_files = input_df[input_filename_column].isin(metadata_df.index)
        if intersect_files.sum() == 0:
            log += "\t\t\t\t###WARNING###: No overlap between samples in metadata and DIA-NN report. Please check whether filename_col is defined correctly.\n\t\t\t\t\t\tReturning not filtered input.\n"
            #sys.exit("###ERROR###: No overlap between samples in metadata and DIA-NN report. please check whether filename_col is defined correctly.")
            return([input_df,log])
        input_df = input_df[intersect_files]
        if len_input_df_before > len(input_df):
            log += "\t\t\t\t###WARNING###: Not all files in input_df are found in metadata. Analysis performed only on " + str(len(intersect_files)) + " files present in metadata.\n"
    return([input_df, log])
'''

def exclude_runs_not_in_metadata(wide_df, metadata_df=None, metadata_filename_col="filename"): 
    """
    Trims input_df (input_df[diann_sample_column]) to same samples present in metadata_df[filename_col].
    author: BZ,OL
    
    Parameters
    ----------
    input_df : pd.DataFrame
        E.g. DIA-NN output report.tsv, long format.
    metadata_df : pd.DataFrame
        Metadata mapping File.Name to shorter run names. The default is None. 
    metadata_filename_col : str
        name of the column in metadata_df specifying filenames (including path) as in DIA-NN report.tsv. The default is "filename".
    
    Returns
    -------
    input_df : pd.DataFrame
        filtered DIA-NN output report.tsv.
    log : str
        string containing messages for report log.
        
    """
    log = ""
    #Exclude runs from input_df not present in metadata
    if metadata_df is not None:
        intersect_files = wide_df.columns.isin(metadata_df[metadata_filename_col])
        if intersect_files.sum() == 0:
            intersect_files = wide_df.columns.isin(metadata_df.index)
        if intersect_files.sum() == 0:
            log += "###ERROR###: No overlap between samples in metadata and DIA-NN report. Please check whether filename_col is defined correctly.\n\t\t\t\t\t\tReturning not filtered input.\n"
            #sys.exit("###ERROR###: No overlap between samples in metadata and DIA-NN report. please check whether filename_col is defined correctly.")
            return wide_df,log
        missing_samples = [sample for sample in wide_df.columns if sample in metadata_df[metadata_filename_col]]
        wide_df = wide_df.drop(columns=missing_samples)
        if len(missing_samples) > 0:
            log += "###WARNING###: Not all files in input_df are found in metadata. Analysis performed only on " + str(len(intersect_files)) + " files present in metadata.\n"
    return wide_df, log


def data_read_column(filename, column_name, dtype=str):
    """
    Reads in one column of a data file (.tsv or .parquet)
    author: OL
    
    Parameters
    ----------
    filename : str
        Path of file to be read in.
    column_name : str
        Column of the file to be read in.
    dtype : dtype
        Data type of the column, e.g. str, int8, float64
        
    Returns
    -------
    column : np.array
        Array containing all the values of the column
    """
    log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tLoading colum "+ column_name +" of input file " + data_df_fp + "\n"
    if "parquet" in filename.split("."):
         column = pd.read_parquet(filename, columns=[column_name], engine="pyarrow")
    else: # read tsv file
         column = pd.read_csv(filename, usecols=[column_name], sep="\t")
    
    try:    
        column = column.astype(dtype).values.reshape(-1)
        log = ""
    except:
        column = column.astype(str).values.reshape(-1)
        log = "###WARNING###: dtype not valid. "+str(column_name)+" was read-in as string.\n"
    return column, log


def select_custom_filtering_columns(ms_filter_columns, ms_filter_operators, ms_filter_values, ms_filter_dtypes, mask_description_df=None, filtering_operation_dict={}):
    """    
    with given custom filters dictionary, ids and values columns this function creates lists with filters for the report.
    author: OL
    
    Parameters
    ----------
    ms_filter_columns : str
        ";"-separated string of filter columns.
    ms_filter_operators : str
        ";"-separated string of filter operators.
    ms_filter_values : str
        ";"-separated string of filter values.
    ms_filter_dtypes : str
        ";"-separated string of filter dtypes.
    mask_description_df: pd.DataFrame, optional
        mapping numbers in mask_df to filters. The default is None, calling creating_filtering_matrix() to create this variable. The default is None.
    filtering_operation_dict : dict, optional
        dictionary containg all filters.
        
    Returns
    -------
    filtering_operation_dict : dict
        dictionary containg all filters.
    log : str
        log for error catching and debugging. 
    """
    
    log = ""
    if mask_description_df is None:
        mask_description_df = pd.DataFrame(columns=["filter","value"])
        
    offset = len(mask_description_df)
    
    filter_columns = str(ms_filter_columns).split(";")
    operators = str(ms_filter_operators).split(";")
    values = str(ms_filter_values).split(";")
    dtypes = str(ms_filter_dtypes).split(";")
    
    if len(filter_columns) != len(operators) or len(filter_columns) != len(values) or len(filter_columns) != len(dtypes):
        log = "###ERROR###: length of filters, operators, values or dtypes not matching.\n"
        return 
    else:
        for ind in range(len(filter_columns)):
            filtering_operation_dict.update({"filter_by_"+filter_columns[ind]: {
                "column": filter_columns[ind].strip(),
                "operator": operators[ind].strip(),
                "value": values[ind].strip(),
                "dtype": dtypes[ind].strip()}
            })
            if filter_columns[ind].strip() not in mask_description_df["filter"].values:
                mask_description_df = pd.concat([mask_description_df,pd.DataFrame(columns=["filter","value"],data = [[filter_columns[ind].strip(),offset+ind]])],axis=0)
    
    return(filtering_operation_dict, mask_description_df, log)



def filtering(filename, filtering_operation_dict, mask_description_df, mask_arr, values_col, keep_track=True):
    """
    Applies filters to values.
    
    Parameters
    ----------
    filename : str
        string specifying DIA-NN output path and filename (report.tsv). The default is None, opening a GUI to select the file.
    filtering_operation_dict : dict
        dictionary containg all filters.
    mask_description_df: pd.DataFrame
        mapping numbers in mask_df to filters.
    mask_arr : np.array
        array in which mask values will be saved
    values_col : str
        values extracted from DIA-NN report.tsv, e.g. Precursor.Quantity.
    keep_track : bool
        if True keeps track of how many entries remain after filtering. The default is True.
    
    Return
    ------
    mask_arr : np.array
        array in which mask values are saved
    log : str
        string containing logs for the step to be written to disk later.
    """
    log = ""
    
    operators = {
    "==": operator.eq,
    "!=": operator.ne,
    ">":  operator.gt,
    ">=": operator.ge,
    "<":  operator.lt,
    "<=": operator.le,
    }
    
    ### NOTE: Overlapping filters (not exclusive)
    if keep_track:
        filtering_tracker_df = pd.DataFrame(index=list(filtering_operation_dict.keys()),columns=["number of excluded precursors (not exclusive)","number of remaining precursors"])    
    else:
        filtering_tracker_df = None
    for key,item in filtering_operation_dict.items():
        try:
            if item["column"] == "__dummy__mass__ms1__":
                filter_arr = data_read_column(filename,item["Ms1.Profile.Corr"],dtype=float)[0]+0.5*data_read_column(filename,item["Mass.Evidence"],dtype=float)[0]
            else:
                filter_arr, log_ = data_read_column(filename,item["column"],dtype=item["dtype"])
                if log_:
                    log += log_
            apply_filter = True
        except:
            if item["column"] == "__dummy__mass__ms1__":
                log += "###WARNING###: Ms1.Profile.Corr and/or Mass.Evidence not in input_df, filter skipped\n"
            else:
                log += "###WARNING###: "+item["column"]+" not in input_df, filter skipped\n"
            apply_filter = False
        if apply_filter:
            try:
                cast_val = np.array(item["value"]).astype(item["dtype"]).item()
                filter_arr = operators[item["operator"]](filter_arr, cast_val)
            except Exception as e:
                log += "###WARNING###: Filtering failed for "+f"{item['column']}{item['operator']}{str(item['value'])}"+", filter skipped.\n"
                log += "###WARNING###: Filtering failed for "+f"{item['column']}{item['operator']}{repr(item['value'])}"+", filter skipped\n"
                log += f"###WARNING###: Reason: {e}\n"
                filter_arr = np.ones(len(filter_arr), dtype=bool)
            
            
            col_ = ""
            
            if item["column"] == "__dummy__mass__ms1__":
                col_ = "Ms1.Profile.Corr + Mass.Evidence"
            elif item["column"] == values_col:
                col_ = "Zero "+values_col
            else:
                col_ = item["column"]
 
            mask_arr[~filter_arr] = mask_description_df.loc[col_].value
            
            if keep_track:
                filtering_tracker_df.loc[key,"number of excluded precursors (not exclusive)"] = sum(mask_arr==mask_description_df.loc[col_].value)
                filtering_tracker_df.loc[key,"number of remaining precursors"] = sum(mask_arr==1.0)
    if keep_track:
        log += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tFiltering completed. "+str(sum(mask_arr!=1.0))+"("+str(np.round(100-sum(mask_arr==1.0)/len(mask_arr)*100,2))+"%) Peptides were excluded during filtering of input_df.\n\n"
        log += filtering_tracker_df.to_string().replace("\n","\n\t") + "\n"
        log += "Note: Filtering per filter is not crisp and can overlap, i.e. precursors can occur multiple times as filtered, thus sum(filtered precursors per step) >= number of total filtered precursors.\n"
    
    return mask_arr, log

def create_mask(mask_arr, mask_description_df, precursor_codes, precursor_labels, file_codes, file_labels, last_filter_value = 14):
    """
    Creates mask and cleans mask description and mask_arr.
    
    Parameters
    ----------
    mask_arr : np.array
        array in which mask values are saved
    mask_description_df: pd.DataFrame
        mapping numbers in mask_df to filters.
    precursor_codes : np.array
        array containg the indices for the rows.
    precursor_labels : np.array
        array containg the labels for the rows matching precursor_codes.
    file_codes : np.array
        array containg the indices for the columns.
    file_labels : np.array
        array containg the labels for the columns matching file_codes.
    last_filter_value : int
        number of last_filter_value in default settings.
    
    Return
    ------
    mask_df : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered
    mask_description_df: pd.DataFrame
        mapping numbers in mask_df to filters. Cleaned for missing entries
    mask_arr : np.array
        Mask array tranfromed to array of bool for elements to keep.
    
    """
    mask_df = pd.DataFrame.sparse.from_spmatrix(coo_matrix((mask_arr, (precursor_codes, file_codes))),
                                                index = precursor_labels,
                                                columns = file_labels).fillna(0.)
    list_of_filters = set(mask_arr)
    value2entry = {mask_description_df.loc[index,"value"]:index for index in mask_description_df.index}
    mask_description_df = mask_description_df.copy()
    for value in range(2,last_filter_value+1):
        if np.float32(value) not in list_of_filters:
            mask_description_df.drop(index = value2entry[np.float32(value)],inplace=True)
    mask_arr = mask_arr == 1
    return mask_df, mask_description_df, mask_arr

def long_to_wide_format(keep_arr, precursor_codes, precursor_labels, file_codes, file_labels, metadata_df, metadata_filename_col="File.Name",metadata_sample_col="Sample.Name", replace_zeroes=True, filter_codes = None):    
    """
    Transforms long format data frame to wide format. 
    author: BZ, OL
    
    Parameters
    ----------
    mask_arr : np.array
        array in which mask values are saved
    precursor_codes : np.array
        array containg the indices for the rows.
    precursor_labels : np.array
        array containg the labels for the rows matching precursor_codes.
    file_codes : np.array
        array containg the indices for the columns.
    file_labels : np.array
        array containg the labels for the columns matching file_codes.
    metadata_df : pd.DataFrame
        metadata mapping diann_sample_column (File.Name/Run) to shorter run names.
    metadata_filename_col : str, optional
        column name with complete filenames in metadata, which match file names in long_df[input_filename_col]. The default is "File.Name".
    metadata_sample_col : str, optional
        column name with (shorter) sample names in metadata. The default is "Sample.Name".
    replace_zeroes : bool, optional
        If True replaces 0s with np.nan. The default is True.
    filter_codes : list of bool, optional
        booloean filter which lines to kepp. The default is None.
        
    Return
    ------
    wide_df : pd.DataFrame
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered and which were measured.
    """
    
    values_arr, log_ = data_read_column(data_df_fp, parameters_dict["input_values_column"], dtype=float)
    if filter_codes is not None:
        values_arr = values_arr[filter_codes]
    values_arr = values_arr.copy()
    values_arr[~keep_arr] = np.nan
    if replace_zeroes:
        values_arr[values_arr==0] = np.nan
    
    wide_df = pd.DataFrame.sparse.from_spmatrix(coo_matrix((values_arr, (precursor_codes, file_codes))),
                                                index = precursor_labels,
                                                columns = file_labels)
    if metadata_filename_col != metadata_sample_col:
        wide_df.columns = metadata_df.reset_index().set_index(metadata_filename_col,drop=False).loc[wide_df.columns,metadata_sample_col]
        wide_df.sort_index(inplace=True,axis=1)
    
    return wide_df

def generate_id_mapping(filename, ids_list, precursor_codes, precursor_suffix_col=None, separator="-", exclude="", log="", filter_codes = None):
    """
    Generates mapping for precursors to proteins to additional data.
    
    Parameters
    ----------
    filename : str
        string specifying MS-dataN output path and filename (e.g. report.tsv)
    ids_list : list of str
        list containing column names to be exctracted-
    precursor_codes : np.array
        array containg the indices for the rows.
    precursor_suffix_col : str
        column_name of precursor suffixes
    separator: str
        separator for suffixes
    exclude : str
        label that should ne be added as a suffix.
    log : str, optional
        string containing the log up to this point.
    filter_codes : list of bool, optional
        booloean filter which lines to kepp. The default is None.
    
    Return
    ------
    id_mapping_df : pd.DataFrame
        dataframe for mapping precursor to protein to other information
    log : str, optional
        string containing the log including warnings raised.
    
    """
    precursor_indices, first_entry = np.unique(precursor_codes, return_index=True)
    most_common_precursor_code = np.bincount(precursor_codes).argmax()
    most_common_precursor_indices = np.where(precursor_codes==most_common_precursor_code)[0]
    
    id_mapping_df = pd.DataFrame(columns=ids_list)
    
    log = ""
    
    if precursor_suffix_col:
        suffixes, _ = data_read_column(filename,precursor_suffix_col,dtype=str)
        if filter_codes is not None:
            suffixes = suffixes[filter_codes]
        mask = (suffixes[first_entry] != "") & (suffixes[first_entry] != exclude)
    else: 
        suffixes = None
    
    for ind,column in enumerate(ids_list):
        try:
            column_, _ = data_read_column(filename,column,dtype=str)
            if filter_codes is not None:
                column_ = column_[filter_codes]
            if ind <= 1 and suffixes is not None:
                id_mapping_df[column] = np.where(mask,np.char.add(np.char.add(column_[first_entry], separator),suffixes[first_entry]),column_[first_entry])
            else:
                id_mapping_df[column] = column_[first_entry]
            if len(set(column_[most_common_precursor_indices])) > 1:
                log += "###WARNING###: ID mapping contained non-unique precursor-column combination for column "+str(column)+". IDs need to be treated with care and additional filters need to be considered.\n"
        except:
            log += "###WARNING###: "+column+" not found in input data. Skipped.\n"
    return id_mapping_df, log
            
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

# DEPRECATED
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
        metadata mapping diann_sample_column (File.Name/Run) to shorter run names.
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
'''    
    
# %%


if __name__ == "__main__":
    
    module = "FILTERING_LONG"
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
    
    #initialize module
    data_df_fp, parameters_dict, id_mapping_df, seed, log, metadata_df, mask_df, mask_description_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = module, 
                                                                                                                                                                  get_metadata=True, 
                                                                                                                                                                  get_id_mapping=False, 
                                                                                                                                                                  get_mask=False, 
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
    
    
    try:
        # Check if loading was successful, else pass error to log and exit   
        helper.logger(log_file,log)
        log = ""
        if metadata_df is None:
            helper.logger(log_file,"###ERROR###: Metadata loading failed.\n")
            if exit_module:
                sys.exit("###ERROR###: Metadata loading failed.")
        if "###WARNING###" in log:
            helper.logger(log_file,"###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"]+"\n")
            if exit_module:
                sys.exit("###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])
        
        # Extract File names and positions in input data
        file_name_column, log_ = data_read_column(data_df_fp, parameters_dict["input_filename_column"], dtype=str)
        file_codes, file_labels = pd.factorize(file_name_column)
        # Clean memory
        del file_name_column
        log += log_
        # Extract Precursor names and positions in input data
        precursor_name_column, log_ = data_read_column(data_df_fp, parameters_dict["precursor_col"], dtype=str)
        if parameters_dict["expand_precursor_name"]:
            precursor_name_column = precursor_name_column.copy()
            precursor_name_suffix, log_ = data_read_column(data_df_fp, parameters_dict["precursor_suffix_col"], dtype=str)
            mask = (precursor_name_suffix != "") & (precursor_name_suffix != parameters_dict["precursor_suffix_exclude"])
            precursor_name_column[mask] = np.char.add(np.char.add([str(p) for p in precursor_name_column[mask]],str(parameters_dict["precursor_suffix_separator"])),[str(s) for s in precursor_name_suffix[mask]])
            del precursor_name_suffix
        precursor_codes, precursor_labels = pd.factorize(precursor_name_column)
        # Clean memory
        del precursor_name_column
        input_id_col = parameters_dict["precursor_col"]
        log += log_
        
        # Rest Mask an Filter array
        mask_arr = np.ones(len(file_codes)).astype("float32") ### FillNa with 0 ### Transform mask to keep array -> trasnform to boolean and filter data ccordingly
        
        # Create mask description
        mask_description_df = pd.DataFrame(columns=["filter","value"])
        mask_description_df["filter"] = ["not measured","measured and kept","Proteotypic","Q.Value","Global.Q.Value",
                                         "GG.Q.Value","PG.Q.Value","Protein.Q.Value","Global.PG.Q.Value","Lib.Q.Value",
                                         "Lib.PG.Q.Value","PEP","Ms1.Profile.Corr","Mass.Evidence","Ms1.Profile.Corr + Mass.Evidence",
                                         "Zero "+parameters_dict["input_values_column"],"Precursor.Charge","sample missingness","precursor missingness",
                                         "Min. "+"XXX"+" Peptides","Min. "+"XXX"+" Precursors","Normalisation Factors","Correction Missingness","Imputed"]
        mask_description_df["value"]=range(0,len(mask_description_df))
        
        # Generate filter operations (column, operator, value, dtype) for DIA-NN
        if parameters_dict["input_origin"]=="diann":
            filtering_operation_dict = select_filtering_columns(filters_dict=parameters_dict, 
                                                                values_col=parameters_dict["input_values_column"],
                                                                input_filename_col=parameters_dict["input_filename_column"])
        else:
            filtering_operation_dict = {}
        
        # Generate filter operations (column, operator, value, dtype) for (additional) custom filters
        if parameters_dict["ms_filter_columns"] is not None:
            filtering_operation_dict, mask_description_df, log_ = select_custom_filtering_columns(ms_filter_columns=parameters_dict["ms_filter_columns"],
                                                                                                  ms_filter_operators=parameters_dict["ms_filter_operators"],
                                                                                                  ms_filter_values=parameters_dict["ms_filter_values"],
                                                                                                  ms_filter_dtypes = parameters_dict["ms_filter_dtypes"],
                                                                                                  mask_description_df = mask_description_df,
                                                                                                  filtering_operation_dict=filtering_operation_dict)
            log += log_
            if "###ERROR###" in log:
                helper.logger(log_file,log)
                if exit_module:
                    sys.exit("###ERROR###: Long format Filtering failed. Please check log file for more details.")
        
        helper.logger(log_file,log)
        
        # Reshape mask_description_df
        mask_description_df.set_index("filter",inplace=True)
        mask_description_df = mask_description_df.astype("float32")
        mask_description_df.loc["Imputed","value"] = 0.99
        mask_description_df.sort_values(by="value",inplace=True)
        
        # Apply filter operations and fill mask array
        mask_arr, log = filtering(data_df_fp, 
                                  filtering_operation_dict, 
                                  mask_description_df, 
                                  mask_arr, 
                                  values_col=parameters_dict["input_values_column"])
        
        # Check if filtering was successful, else pass error to log and exit
        if sum(mask_arr==1) > 0 and "###ERROR###" not in log:
                log += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tLong format filtering successful.\n"
                helper.logger(log_file,log)
                log=""
        elif sum(mask_arr==1) == 0:
            log += "\t\t\t\t###WARNING###: Filtered input_df is empty. Please check if columns on which filters are based are not empty (sometimes depending on library used e.g. Genes and therefore Q-Values on gene-level are missing), else adjust filtering parameters.\n"
            helper.logger(log_file,log)
            if exit_module:
                sys.exit("###ERROR###: Long format Filtering failed.")
        elif "###ERROR###" in log:
            log += "\n###ERROR###: Long format Filtering failed."
            helper.logger(log_file,log)
            if exit_module:
                sys.exit("###ERROR###: Long format Filtering failed. Please check the log for more details.")
        
        # Remove missing entries
        if any(precursor_codes == -1):
            filter_codes = precursor_codes != -1
            file_codes = file_codes[filter_codes]
            precursor_codes = precursor_codes[filter_codes]
            mask_arr = mask_arr[filter_codes]
        else:
            filter_codes = None
        
        # Generate mask, clean mas description and tranform mask_arr to boolean array
        mask_df, mask_description_df, mask_arr = create_mask(mask_arr, mask_description_df, 
                                                             precursor_codes, precursor_labels, 
                                                             file_codes, file_labels, 
                                                             last_filter_value = 14)
        
        # Transform long format to wide format
        wide_df = long_to_wide_format(mask_arr, 
                                      precursor_codes, precursor_labels, 
                                      file_codes, file_labels, 
                                      metadata_df, 
                                      metadata_filename_col=parameters_dict["metadata_filename_column"],
                                      metadata_sample_col=parameters_dict["metadata_sample_column"],
                                      filter_codes = filter_codes)
        
        # Clean memory
        del mask_arr
        del file_codes
        del file_labels
        
        # Extract variables from parameters, which need adjustment
        ids_list = [parameters_dict["precursor_col"],parameters_dict["protein_col"]]
        
        # Get columns to be saved in ID-mapping
        if parameters_dict["input_origin"]=="diann":
            if parameters_dict["input_ids_columns"] is not None:
                ids_list = ids_list + [col for col in ["Stripped.Sequence","Genes","Proteotypic"] if col not in [parameters_dict["precursor_col"],parameters_dict["protein_col"]]] + [col for col in parameters_dict["input_ids_columns"].split(";") if col not in ["Stripped.Sequence","Genes","Proteotypic"]+[parameters_dict["precursor_col"],parameters_dict["protein_col"]]]
            else:
                ids_list = ids_list + [col for col in ["Stripped.Sequence","Genes","Proteotypic"] if col not in [parameters_dict["precursor_col"],parameters_dict["protein_col"]]]
        else:
            if parameters_dict["input_ids_columns"] is not None:
                ids_list = ids_list + [col for col in parameters_dict["input_ids_columns"].split(";") if col not in [parameters_dict["precursor_col"],parameters_dict["protein_col"]]]
        
        # Generate ID-mapping
        if parameters_dict["expand_precursor_name"]:
            id_mapping_df, log = generate_id_mapping(data_df_fp, ids_list, precursor_codes, parameters_dict["precursor_suffix_col"], parameters_dict["precursor_suffix_separator"], parameters_dict["precursor_suffix_exclude"], log = log, filter_codes = filter_codes) 

        else:
            id_mapping_df, log = generate_id_mapping(data_df_fp, ids_list, precursor_codes, log = log, filter_codes = filter_codes) 
        
        # Clean memory
        del precursor_codes
        del precursor_labels
        
        # Add Peptide column
        if parameters_dict["peptide_col"] not in id_mapping_df.columns:
            id_mapping_df[parameters_dict["peptide_col"]] = create_peptide_ids(id_mapping_df[parameters_dict["precursor_col"]])
        
        # Rename missing protein names with other column
        ### Not Tested yet
        #if pd.isna(id_mapping_df[parameters_dict["protein_col"]]).sum() > 0:
        #    id_mapping_df = name_unknown_proteins(id_mapping_df=id_mapping_df,
        #                                          input_fp=data_df_fp, 
        #                                          ids_list=ids_list, 
        #                                          protein_col=parameters_dict["protein_col"], 
        #                                          replace_id_col=parameters_dict["alternative_protein_col"])    
            
        log += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tID mapping successfully generated.\n"
        helper.logger(log_file,log)
        
        # Sample filtering (trimming to samples specified by metadata_df)
        wide_df, log = exclude_runs_not_in_metadata(wide_df=wide_df, 
                                                     metadata_df=metadata_df, 
                                                     metadata_filename_col=parameters_dict["metadata_sample_column"])
        wide_df = wide_df.sparse.to_dense()
        wide_df.index.names = [parameters_dict["precursor_col"]]
        wide_df.dropna(how="all", inplace=True)
        helper.logger(log_file,log)

        # Rename filter_mask columns
        mask_df = helper.rename_filter_mask_columns(mask_df, 
                                                    metadata_df, 
                                                    filename_col=parameters_dict["metadata_filename_column"], 
                                                    sample_col=parameters_dict["metadata_sample_column"])
        
        mask_df = mask_df.sparse.to_dense()
        mask_df.index.names = [parameters_dict["precursor_col"]]
        
        # Log Transformation (Base 2)
        wide_df = helper.log_transform(wide_df, base=2)
        log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tWide format data logtransformed to base 2\n"
        #helper.logger(log_file,log)
        
        # Save data
        id_mapping_df.to_csv(os.path.join(out_dir_fp, "id_mapping.tsv"), index=False, sep="\t")
        log += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tID-Mapping saved successfully to " + os.path.join(out_dir_fp, "id_mapping.tsv")+"\n"
        helper.close_module(data_processed_df=wide_df, masked_df=mask_df, mask_description_df=mask_description_df, log_file=log_file, module_output=log, out_dir_fp=out_dir_fp, outname=outname, mask=True, exit_module=exit_module)
        
    except Exception as e:
             e_type, e_obj, e_tb = sys.exc_info()
             script_name = os.path.split(e_tb.tb_frame.f_code.co_filename)[1]
             helper.logger(log_file, "###ERROR### Module " +script_name + " ran into an error at line : " + str(e_tb.tb_lineno) +"\n\t\t\t\t" + str(e) + "\n")
             if exit_module:
                 sys.exit("Filtering_long exited with an error.")
