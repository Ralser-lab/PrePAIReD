"""
Input: X
Output: X

functions:
	- zscore_transformation: calculates z-scores based on mean and sd or median and median absolute deviation (MAD)
	- log_transform: applies a logarithmic transformation to a DataFrame with the specified base
	- unlog_transform: transforms intensities of wide-format precursor or protein intensities over samples from log-space to not-log-space for a specified base
	- mad: calculate median absolute deviation
	- (calculate_mean_over_refs: calculate a mean sample over reference samples specified via refstring and exclude original samples specified via refstring)
	- (calculate_log2FC_to_ref: calculates log2 fold changes to a reference sample specified by ref)
	- (calculate_log2foldchange: calculate the log2 fold change for each protein/precursor relative to the median (or mean) intensity across all samples)
	- class OverwriteException: class for raising OverwriteException
	- create_directory: creates (if not existent) a directory
	- replace_invalid_characters: replaces characters that cannot be decoded in standard encoding, e.g. replace greek letters
	- create_logfile: creates a logfile following changes during pre-processing => DO WE NEEDED?
	- create_filtering_matrix: creates a pandas.DataFrame keeping track for every precursor in every run why it was filtered, imputed, ...
	- rename_filter_mask_columns: renames diann_sample_column ("File.Name" or "Run") in filter mask columns to shorter sample names
	- create_filename_annotations_file: automatically create metadata with mappings between diann_sample_column, shorter run names and sample group
	- (check_filename_annotations_on_unique_indices: control if metadata has unique entries)
	- load_metadata: load the metadata file (tsv or Excel)
	- weighted_median: calculate the weighted median for an array with weights W
	- weighted_nanmedian: calculate weighted median ignoring NaNs and zero weights
	- (insert_column_level: inserts a new level to MultiIndex columns)
	- (annotate_metadata_to_profiles: adds all metadata as multi-layer column names to profiles)
	- create_imputation_table: create a template for an involved table to define different impuation methods across sample and precursor hierarchies
	- logger: iteratively callable function to add information to the existing logfile of the current preprocessing run
	- get_parser: generates the parser to parse arguments from the terminal to the Python Code
	- get_parser_arguments: extraction of parser arguments
	- load_data: loads data for processing
	- load_parameters: read in parameter *.json file
	- load_id_mapping: reads in ID mapping
	- get_seed: generate the seed for the module
	- add_parameters_to_log: writes parameters for the module to the log file
	- load_mask: loads the filtering/imputation mask
	- load_mask_description: loads description for the mask identifers
	- save_mask: saves (updated) mask to disk
	- save_data: saves data after processing to disk
	- init_module: initializes the module and reads in all necessary data and arguments
	- close_module: saves data to disk and exits the module

"""

import re
import numpy as np
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import os.path
import os
import glob
import time
import argparse
import json
import sys
import itertools

def zscore_transform(df, robust=False, axis:int=0):
    """

    Calculates z-scores based on mean and sd or median and median absolute 
    deviation (MAD) on the preferred axis (0:rows, 1: columns).

    if robust == True use: z-scores over samples as (X-median(X))/MAD(X)
    
    author: BZ, SB, CD
    
    Parameters
    ----------
    data : pd.DataFrame
        dataframe with features as columns and samples as rows
    axis : int, optional
        Set to 0 if columns in dataframe are features and samples are rows, set 
        to 1 if rows and columns are transposed. The default is 0.

    Returns
    -------
    zdata : pd.DataFrame
        dataframe with same shape as data and z-scores as values

    """

    # Throw a ValueError if axis is not 1 or 0
    if axis not in [0, 1]:
        raise ValueError("Axis must be either 0 (rows) or 1 (columns)")

    if robust:
        mads=df.sub(df.median(axis=axis),axis=1-axis).abs().median(axis=axis)
        df_z=df.sub(df.median(axis=axis),axis=1-axis).div(mads,axis=1-axis)
    else:
        means = df.mean(axis=axis)
        sds = df.std(axis=axis)
        df_z = df.sub(means,axis=1-axis).div(sds, axis=1-axis)

    return df_z 

def log_transform(
        df, base:int=2, 
        add_constant=False, constant=1):
    
    """
    Applies a logarithmic transformation to a DataFrame with the specified base.

    Authors: SB, BZ, CD

    Parameters:
    -----------
    df : pandas.DataFrame
        The input DataFrame containing numerical values to be transformed.
        
    base : int (default=2)
        The base of the logarithm. 
        
    add_constant : bool, optional (default=False)
        If True, a constant value will be added to all elements of the DataFrame 
        before applying the logarithm transformation. This can be useful if your 
        data contains zero or negative values (since log of zero or negative is 
        undefined).
        
    constant : int, float, optional (default=1)
        The constant value to be added to each element of the DataFrame, used 
        only if `add_constant=True`.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame with the log-transformed values. The transformed values 
        will have the same shape as the input `df`.
        Missing values (NaN) will remain NaN in the output DataFrame.

    Raises:
    -------
    ValueError:
        If the DataFrame contains non-numeric values or if the base is 
        non-positive or equal to 1.
    """

    # Check for invalid base values
    if base <= 0 or base == 1:
        raise ValueError("Base must be greater than 0 and not equal to 1")
    # Get NaNs
    nan_mask = df.isna()
    # Add constant
    if add_constant:
        df = df + constant
    # Log transform
    df = np.log(df) / np.log(base)
    # Replace Nan's
    df[nan_mask] = np.nan
    return df


def unlog_transform(
        df, base=2):

    """
    Transforms intensities of wide-format precursor or protein intensities 
    over samples from log_{Base}-space to not-log_{Base}-space.
        
    Author: BZ
    
    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame where rows could represent precursor or protein IDs and 
        columns represent samples. 

    base : int, optional (default=2)
        The base of the logarithmic transformation. By default, it is set to 2 
        for log2-space. You can adjust this to match the base of the original log 
        transformation (e.g., 10 for log10).
    
    Returns
    -------
    df: pd.DataFrame
        A DataFrame of the same shape as the input, but with values transformed 
        back to the original scale by applying the exponentiation with the 
        specified base.

    """
	# Check for invalid base values
    if base <= 0 or base == 1:
        raise ValueError("Base must be greater than 0 and not equal to 1")
    df = base**df
	# Replace 1 with 0 to preserve zeros before transformation
    df = df.replace(1,0)
    return(df)

def mad(x):
    """
    Calculate median absolute deviation of a vector x
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        Array with float values.

    Returns
    -------
    float 

    """
    return np.nanmedian(np.abs(x-np.nanmedian(x)))

"""
BZ: 
functions used for ratio-based analysis mainly (still normalization functions)
e.g. normalization of sample fold changes, normalization to a mean
might be excluded for now
"""

def calculate_mean_over_refs(
        df, refstring=r'^BY4741.*'):
    """
    Calculate a mean sample over reference samples specified via refstring and 
    exclude original samples specified via refstring.
    
    author: BZ
    
    Parameters
    ----------
    df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.
    refstring : str, optional
        regular expression used to identify reference samples. 
        The default is r'^BY4741.*'.

    Returns
    -------
    df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.

    """
    #define reference columns
    r=re.compile(refstring)
    refs=[m.group(0) for m in map(r.match,df.columns) if m is not None]
    #calculate mean over reference
    df["ref"]=df[refs].mean(axis=1)
    #remove ref strains
    df.drop(refs,axis=1,inplace=True)
    return df


def calculate_log2FC_to_ref(df,ref ,level=None):
    """
    Calculates log2 fold changes to a reference sample specified by ref.
    
    author: BZ
    
    Parameters
    ----------
    df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions.
    ref : str
        column name of reference sample to which fold changes should be 
        calculated.
    level : str, optional
        For multi-column df please provide a level in which ref is found. 
        The default is None.

    Returns
    -------
    log2fc_df : pandas.DataFrame
        rows equal proteins/precursors, columns equal samples/runs/fractions, 
        values are log2 fold changes to reference column.

    """
    if level is None:
        log2fc_df = np.log2(df.div(df[ref], axis=0))
    else:
        log2fc_df = np.log2(
            df.div(df.xs(ref, axis=1,level=level).iloc[:,0].tolist(), axis=0))
    return log2fc_df 


def calculate_log2foldchange(df, robust=True):
    """
    Calculate the log2 fold change for each protein/precursor relative to the 
    median (or mean) intensity across all samples.
    
    author: BZ
    
    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame where rows represent proteins/precursors, and columns 
        represent samples/runs/fractions. The values in the DataFrame are 
        intensity measurements.

    robust : bool, optional (default=True)
        If True, the log2 fold change is calculated relative to the median 
        intensity across samples for each protein. If False, the fold change 
        is calculated relative to the mean intensity across samples.

    Returns
    -------
    pd.DataFrame
        A DataFrame with the same shape as the input, where each value is the 
        log2 fold change of the corresponding protein/precursor relative to 
        the median (or mean) intensity across all samples. Zero values are 
        replaced with `NaN` during calculation and restored to 0 in the result.

    """
    df.replace(0, np.nan, inplace=True)
    if robust:
        log2fcs = np.log2(df.div(df.median(axis=1),axis=0))
    else:
        log2fcs = np.log2(df.div(df.mean(axis=1),axis=0))
    
    log2fcs = log2fcs.replace(np.nan,0)
    return log2fcs 
    


class OverwriteException(Exception):
	def __init__(self,message):
		super().__init__(message)


'''
in inheriting file

from OverwriteException import OverwriteException

try:
	if overwrite==True:
		raise OverwriteException("overwriting of existing file requested")
	load("try") # not reached if OverwriteException
except (FileNotFoundError,OverwriteException):
	pass

'''




def create_directory(directory_name):
    """
    creates (if not existent) the directory for analysis
    
    author: BZ
    
    Parameters
    ----------
    directory_name : str
        name and path of directory.

    Returns
    -------
    None.

    """
    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)
        print("folder "+directory_name+" created")
    return



def replace_invalid_characters(filename,sep="\t",to_be_replaced="µ",replace_with="u",encoding="cp1253",appendix=""):
    """
    Replaces characters that cannot be decoded in standard encoding, e.g. replace greek letters.
    
    author: BZ
    
    Parameters
    ----------
    filename : str
        Name of the dataframe to read.
    sep: str, optional
        separator passed to pandas.read_csv as argument sep. The default is "\t"
    to_be_replaced : str, optional
        pattern that needs to be replaced. The default is "µ".
    replace_with : str, optional
        pattern to replace to_be_replaced with. The default is "u".
    encoding : str, optional
        encoding for reading dataframe with invalid characters. The default is "cp1253".
    appendix : str, optional
        appendix to add to filename before saving. The default is "", overwriting the initial filename with the corrected version.
    Returns
    -------
    None.

    """
    df=pd.read_csv(filename,sep=sep,encoding=encoding)
    df=df.replace(to_be_replaced,replace_with,regex=True)
    df.to_csv(filename[:-4] + appendix + filename[-4:],sep=sep,index=False)
    return


#saves parameters and progress to log file
def create_logfile(logfilename,ids,savefilename,analysis_folder,profile_annotation_file,values,filters,filename):
    """
    Creates a logfile following changes during pre-processing.
    
    author: BZ
    
    Parameters
    ----------
    logfilename : str
        name (path + filename) of the logfile
    ids : list
        list of ids extracted from DIA-NN report.tsv.
    savefilename : str
        name of pre-processed file.
    analysis_folder : str
        name of directory where analysis (pre-processed file), logfile, metadata are stored.
    profile_annotation_file : str
        name of metadata file.
    values : str
        values extracted from DIA-NN report.tsv, e.g. Precursor.Quantity.
    filters : dict
        dictionary containing the used filters.
    filename : str
        name of DIA-NN report.tsv.

    Returns
    -------
    None.

    """
    #create logfile
    logfile=open(logfilename,"w")
    #write some general information
    logfile.write("Logfile for creation of "+savefilename+".tsv\nResults are being saved to "+analysis_folder+"\nThe filtering also generates following files:\n\t- "+profile_annotation_file + " with meta data annotations to the file (shorter sample names and sample group for every sample\n\t- "+logfilename.split("/")[-1]+" (this file)\n\t- "+savefilename+"_filtermask.tsv"+" keeping track why precursors were excluded or imputed\n\t- "+savefilename+"_filtermask_description.tsv"+" providing a mapping between the values in filter_mask and respective filters\n")
    logfile.write("Analysis on DIA-NN output: "+filename+"\n")
    logfile.write("\n-----------------------------------------------\nPARAMETERS\n\n")
    logfile.write("ID columns in "+savefilename+":\t\t"+", ".join(ids)+"\n")
    logfile.write("Values displayed in "+savefilename+":\t"+values+"\n")
    logfile.write("\n-----------------------------------------------\n")
    logfile.write("\nFiltering Parameters\n")
    #write all filters and filtering parameters
    for key,value in filters.items():
        if len(key)<15:
            logfile.write(key+":\t\t\t\t"+str(value)+"\n")
        elif len(key)<22:
            logfile.write(key+":\t\t\t"+str(value)+"\n")
        elif len(key)<28:
            logfile.write(key+":\t\t"+str(value)+"\n")
        else:
            logfile.write(key+":\t"+str(value)+"\n")
    if filters["filter_by_precursorcharge"]:
        logfile.write("Minimum precursor charge:\t\t"+str(filters["preccharge_min"])+"\nMaximum precursor charge:\t\t"+str(filters["preccharge_max"])+"\n")
    logfile.write("\n-----------------------------------------------\n")
    logfile.close()
    return
    

#create a matrix with precursors as rows and files/runs/samples/fractions as columns to keep track of filtered and imputed precursors
def create_filtering_matrix(file,metadata_df,min_peptides,min_precursors,values,ids="Precursor.Id",filename_col="filename",diann_sample_column="File.Name"):
    """
    Creates a pandas.DataFrame keeping track for every precursor in every run why it was filtered, imputed, ...
    
    author: BZ
    
    Parameters
    ----------
    file : pd.DataFrame
        DIA-NN output report.tsv.
    metadata_df : pd.DataFrame
        Metadata mapping File.Name to shorter run names.
    min_peptides : int
        minimum number of peptides required to quantify a protein.
    min_precursors : int
        minimum number of precursors required to quantify a protein.
    values : str
        values extracted from DIA-NN report.tsv, e.g. Precursor.Quantity.
    ids : list, optional
        list of ids extracted from DIA-NN report.tsv. The default is ["Precursor.Id","Stripped.Sequence","Protein.Group","Genes"].
    filename_col : str, optional
        column name with complete filenames in metadata. The default is "filename".
    diann_sample_column : str, optional
        name of the sample names in the DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
        
    Returns
    -------
    mask_df : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered
    mask_description_df: pd.DataFrame.
        mapping numbers in mask to filters

    """
    log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tCreating filtering matrix\n"
    #get sorted File.Names
    #samples=[x.replace("\\","/").split("/")[-1] for x in file[diann_sample_column].unique()]
    metadata_df.set_index(filename_col,inplace=True,drop=False)
    samples=[x for x in file[diann_sample_column].unique()]
    if len(set(samples)) > len(set(metadata_df.index)):
        log += "\t\t\t###WARNING###: Not all samples in DIA-NN output annotated in metadata\n"
    elif len(set(samples)) < len(set(metadata_df.index)):
        log += "\t\t\t###WARNING###: Not all samples in metadata found in DIA-NN output\n"
    #get intersection between filenames and metadata
    samples=list(set(samples) & set(metadata_df.index))
    #sort File.Names
    samples.sort()
    columns_filenames=metadata_df.loc[samples,filename_col]
    mask_df = pd.DataFrame(index=file.set_index(ids).index.unique(), columns=columns_filenames,data=0)
    # dataframe explaining reasons why peptide was excluded
    mask_description_df = pd.DataFrame(columns=["filter","value"])
    #mask_description["filter"] = ["not measured","measured and kept", "excessing_NaNs","Q.Value","Global.Q.Value","GG.Q.Value","PG.Q.Value","Protein.Q.Value","Global.PG.Q.Value","Lib.Q.Value","Lib.PG.Q.Value","Precursor.Charge","Proteotypic","Min. "+str(min_peptides)+" Peptides","Zero "+values,"Imputed"]
    #reordered to same order as filtering steps
    mask_description_df["filter"] = ["not measured","measured and kept","Proteotypic","Q.Value","Global.Q.Value","GG.Q.Value","PG.Q.Value","Protein.Q.Value","Global.PG.Q.Value","Lib.Q.Value","Lib.PG.Q.Value","PEP","Ms1.Profile.Corr","Mass.Evidence","Ms1.Profile.Corr + Mass.Evidence","Zero "+values,"Precursor.Charge","sample missingness","precursor missingness","Min. "+str(min_peptides)+" Peptides","Min. "+str(min_precursors)+" Precursors","Imputed"]
    mask_description_df["value"]=range(0,len(mask_description_df))
    #write mask description into logfile
    log += "\nFiltering matrix to keep track over excluded precursors and proteins created.\nReasons for exclusion of precursors as specified in filter mask\n" + mask_description_df.to_string(index=False) +"\n-----------------------------------------------\n\n"
    mask_description_df.set_index("filter",inplace=True)
    mask_description_df.loc["Imputed","value"] = 0.99
    return([mask_df,mask_description_df,log])



def rename_filter_mask_columns(filter_mask,metadata_df,filename_col="filename",sample_col="sample"):
    """
    Renames diann_sample_column ("File.Name" or "Run") in filter mask columns to shorter names
    
    author: BZ
    
    Parameters
    ----------
    filter_mask : pd.DataFrame
        output filtering matrix from create_filtering_matrix().
    metadata_df : pd.DataFrame
        Metadata mapping DIA-NN samples to shorter run names.
    filename_col : str, optional
        column name with complete filenames in metadata. The default is "filename".
    sample_col : str, optional
        column name with sample names in metadata. The default is "name on x-axis"

    Returns
    -------
    filter_mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.

    """
    summarized=None
    if "Precursor summarized for Protein" in filter_mask.columns:
        summarized=filter_mask["Precursor summarized for Protein"]
        filter_mask.drop("Precursor summarized for Protein",axis=1,inplace=True)
    #trim filter_mask on samples present in metadata
    samples=list(set(metadata_df[filename_col]) & set(filter_mask.columns))
    filter_mask=filter_mask[samples]
    filter_mask.columns=metadata_df.set_index(filename_col,drop=False).loc[filter_mask.columns,sample_col]
    #filter_mask.columns.names=["sample","group"] ### added 
    #sorting columns
    filter_mask.sort_index(axis=1,level=[1,0],inplace=True)
    if summarized is not None:
        filter_mask["Precursor summarized for Protein"]=summarized
    return(filter_mask)





# Annotation of DIA-NN output File.Names
# profile_order specifies order of quant files in output
def create_filename_annotations_file(path, file, profile_annotation_file, profile_order=None, raw_data_folder=None, profile_annotation="", filename2run_pattern=None, filename2samplegroup_pattern="1", overwrite_filename_annotations=False, filename_col="filename", sample_col="name on x-axis", sample_group_col="sample group", encoding="utf-8", diann_sample_column="File.Name", logfilename=None):
    """
    automatically create metadata with mappings between diann_sample_column, shorter run names and sample group.
    Recommendation: set profile_annotation to "", provide if possible filename2run_pattern and filename2samplegroup_pattern to get shorter sample and group names.
    If neccessary, edit output with MS Excel.
    
    author: BZ
    
    Parameters
    ----------
    path : str
        path of DIA-NN output file.
    file : pd.DataFrame
        (filtered) DIA-NN output report.tsv.
    profile_annotation_file : str
        name of metadata file.
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
    filename_col : str, optional
        column name with complete filenames in metadata. The default is "filename".
    sample_col : str, optional
        column name with sample names in metadata. The default is "name on x-axis"
    sample_group_col: str, optional
        column name with sample group names in metadata. The default is "sample group". 
    encoding: str, optional
        the encoding passed to pandas.read_csv. The default is "utf-8".
    diann_sample_column : str, optional
        name of the sample names in the DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    logfilename : str
        name (path + filename) of the logfile. The default is None, disabling logging.
        
    Raises
    ------
    OverwriteException
        Overwriting of file requested by user (overwrite_filename_annotations=True).

    Returns
    -------
    metadata_df : pd.DataFrame
        automatically created metadata.

    """
    # if order of files is given, sorts files from raw_data_folder according to this order, if none returns [] and continues
    if profile_order is None:
        if raw_data_folder is None:
            raw_data_folder = path+"raw_data/"
        # profile_order=file["File.Name"].unique().tolist() #only with occurences
        files_list = glob.glob(path+raw_data_folder+"*.quant")
        # only take filename without path and without filetype ending
        profile_order = list(map(lambda x: x.split("\\")[-1][:-6], files_list))
        profile_order.sort()
    try:  # try to load existing file profile_annotation_file or if requested open GUI to select file
        if profile_annotation == "from_file":
            root = Tk()
            root.withdraw()
            profile_annotation_file_GUI = askopenfilename(initialdir="../")
            root.destroy()
            metadata_df = pd.read_csv(profile_annotation_file_GUI, sep="\t", dtype=str, encoding=encoding)
        else:
            if overwrite_filename_annotations == False:
                metadata_df = pd.read_csv(profile_annotation_file, sep="\t", dtype=str, encoding=encoding)
                print("Annotations to raw data filenames found and loaded")
                if logfilename is not None:
                    logfile=open(logfilename,"a")
                    logfile.write("\nAnnotations to raw data filenames found and loaded\n")
                    logfile.close()
            else:
                raise OverwriteException("overwriting of existing file requested")
                if logfilename is not None:
                    logfile=open(logfilename,"a")
                    logfile.write("\nOverwriting of annotations to raw data filenames requested by user.\n")
                    logfile.close()
    except (FileNotFoundError, OverwriteException):
        # create dataframe with mappings between original filenames and shorter versions for further analysis
        metadata_df = pd.DataFrame(
            columns=[filename_col, sample_col, sample_group_col])
        if profile_order != []:
            metadata_df[filename_col] = profile_order
        else:
            metadata_df[filename_col] = file[diann_sample_column].unique()
        if profile_annotation == "manual":
            print("Please enter name for every file")
            for i in metadata_df[filename_col]:
                print("filename: "+str(i))
                name = input()
                # filename != index
                metadata_df.loc[i, sample_col] = name
        elif profile_annotation == "Excel":
            metadata_df.to_csv(
                profile_annotation_file, index=False, sep="\t")
        elif profile_annotation == "IEF":
            print("Please enter first pH-value (only enter numbers)")
            start = float(input())
            print("Please enter last pH-value (only enter numbers)")
            end = float(input())
            metadata_df[sample_col] = np.linspace(
                start, end, len(profile_order))
        elif profile_annotation == "SEC":
            print("Please enter first size (only enter numbers)")
            start = float(input())
            print("Please enter last size (only enter numbers)")
            end = float(input())
        elif profile_annotation == "numbers":
            metadata_df[sample_col] = range(
                1, len(profile_order)+1)
        else:  # no names specified, take filenames
            if filename2run_pattern is None:
                metadata_df[sample_col] = metadata_df[filename_col]
            else:
                metadata_df[sample_col] = [eval(filename2run_pattern) for x in metadata_df[filename_col]]
        metadata_df[sample_col] = metadata_df[sample_col].astype(
            str)
        # save profile annotation file
        if profile_annotation_file is not None:
            metadata_df[sample_group_col] = [eval(filename2samplegroup_pattern) for x in metadata_df[filename_col]]
            metadata_df.to_csv(
                profile_annotation_file, index=False, sep="\t")
        print("Annotations to raw data created and saved")
        if logfilename is not None:
            logfile=open(logfilename,"a")
            logfile.write("Annotations to raw data created and saved")
            logfile.close()

    if "FileName_trimmed" not in metadata_df.columns:
        metadata_df["FileName_trimmed"] = list(map(lambda x: x.replace(
            "/", "\\").split("\\")[-1], metadata_df[filename_col].tolist()))
    metadata_df.set_index("FileName_trimmed", inplace=True)
    if profile_annotation_file is not None and overwrite_filename_annotations:
        metadata_df.to_csv(
           profile_annotation_file, index=True, sep="\t")
        print("Annotations ro raw data overwritten")
    return(metadata_df)



def check_filename_annotations_on_unique_indices(profile_annotation_file, sample_col="name on x-axis", sample_group_col="sample group"):
    """
    Control if metadata has unique entries
    
    author: BZ
    
    Parameters
    ----------
    profile_annotation_file : str
        name of metadata file.
    sample_col : str, optional
        column name with sample names in metadata. The default is "name on x-axis"
    sample_group_col: str, optional
        column name with sample group names in metadata. The default is "sample group". 
        
    Returns
    -------
    None.

    """
    metadata_df = pd.read_csv(profile_annotation_file, sep="\t", dtype=str)
    if len(set(metadata_df[[sample_col,sample_group_col]].itertuples(index=False,name=None))) < len(metadata_df):
        print(metadata_df[metadata_df[[sample_col,sample_group_col]].duplicated()][[sample_col,sample_group_col]])
        input("Correct duplicated indices and press Enter...")
    return



def load_metadata(metadata_fp,metadata_sheet="Sheet1", filename_col="File.Name", encoding="utf-8"):
    """
    load the metadata file
    
    author: BZ
    
    Parameters
    ----------
    metadata_fp : str
        name of the metadata file. (The default is None, if None a GUI opens to select the file.)
    metadata_sheet : str, optional
        if metadata is an excel file, metadata_sheet specifies the sheetname. The default is "Sheet1".
    filename_col : str, optional
        column name with complete filenames in metadata. The default is "File.Name".
    encoding : str, optional
        encoding passed to pandas read_csv

        
    Returns
    -------
    metadata_df : pd.DataFrame
        the metadata file as dataframe
    log : str
        string for report log messages

    """
    '''
    if metadata_fp is None:
        print("Please select file containing meta data.")
        root = Tk()
        root.withdraw()
        metadata_fp = askopenfilename(
            initialdir="S:/AG/AG-CF-HTMS/CoreFacility/")
        root.destroy()
    '''
    log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tLoading metadata " + metadata_fp + "\n"
    # read tsv file
    if metadata_fp[-4:] in [".txt",".tsv"]:
        metadata_df = pd.read_csv(metadata_fp, sep="\t", encoding=encoding)
    elif metadata_fp[-4:]==".csv":
        metadata_df = pd.read_csv(metadata_fp, sep=",", encoding=encoding)
    #read xlsx file
    elif metadata_fp[-5:]==".xlsx":
        metadata_excel=pd.ExcelFile(metadata_fp)
        #get sheets, check if sheet in there
        if metadata_sheet not in metadata_excel.sheet_names:
            log += "###ERROR### Metadata file " + metadata_fp + " does not contain provided metadata_sheet " + metadata_sheet
            return([None,log])
        metadata_df=metadata_excel.parse(metadata_sheet)
        metadata_excel.close()
    else:
        log += "\t\t\t###ERROR### Please implement how how to read metadata of filetype "+metadata_fp[-5:] + "\n"
        return([None,log])
    #exclude all annotations without filename (e.g. not measured yet)
    if filename_col in metadata_df.columns:
        metadata_df=metadata_df.dropna(subset=[filename_col])
    else:
        log += "\t\t\t###WARNING### Column " + filename_col +" not found in metadata columns."
    #if "FileName_trimmed" not in metadata_df.columns and filename_col in metadata_df.columns:
    #    metadata_df["FileName_trimmed"] = list(map(lambda x: x.replace("/", "\\").split("\\")[-1], metadata_df[filename_col].tolist()))
    #if "FileName_trimmed" in metadata_df.columns:
    #    metadata_df.set_index("FileName_trimmed", inplace=True)
    # drop duplicated indices
    if sum(metadata_df.index.duplicated())>0:
        log += "\t\t\t###WARNING### " + str(sum(metadata_df.index.duplicated()))+" indices duplicated in metadata. Duplicated indices are being dropped.\n"
        metadata_df=metadata_df.drop_duplicates()
    log += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tMetadata loaded successfully.\n"
    return([metadata_df, log])




#source: https://www.geeksforgeeks.org/program-to-find-weighted-median-of-a-given-array/ , access: 02.11.2022 12:00
def weighted_median(arr,W):
    """
    Calculate the weighted median for an array with weights W.
    author: BZ
    
    Parameters
    ----------
    arr : np.array
        array for which the median should be calculated.
    W : np.array
        weights for every element of arr.

    Returns
    -------
    weighted median of the array arr with weights W.

    """
    # Store pairs of arr[i] and W[i]
    pairs = []
    for index in range(len(arr)):
        pairs.append([arr[index], W[index]])
    # Sort the list of pairs w.r.t.
    # to their arr[] values
    pairs.sort(key = lambda p: p[0])
    # If N is odd
    if len(arr) % 2 != 0:
        # Traverse the set pairs
        # from left to right
        sums = 0
        for element, weight in pairs:
            # Update sums
            sums += weight
            # If sum becomes > 0.5
            if sums > 0.5:
                return(element)
    # If N is even
    else:
        # For lower median traverse
        # the set pairs from left
        sums = 0
        for element, weight in pairs:
            # Update sums
            sums += weight
            # When sum >= 0.5
            if sums >= 0.5:
                lower_median=element
                break
        # For upper median traverse
        # the set pairs from right
        sums = 0
        for index in range(len(pairs)-1, -1, -1):
            element = pairs[index][0]
            weight = pairs[index][1]
            # Update sums
            sums += weight
            # When sum >= 0.5
            if sums >= 0.5:
                upper_median=element
                break
        return(np.mean([lower_median,upper_median]))
    


def weighted_nanmedian(values,weights):
    """
    Calculate weighted median ignoring NaNs and zero weights, scales weights to a sum of 1, returns nan for weighted mean of empty arrays
    author: BZ
    
    Parameters
    ----------
    values : list or np.array
        list or array for which the weighted median should be calculated.
    weights : list or array
        list or array providing the weights for values.

    Returns
    -------
    weighted median of the array arr with weights W ignoring missing values.

    """
    #remove missing values and convert input to np.arrays
    nans=np.isnan(values)
    values=np.asarray(values)[~nans]
    weights=np.asarray(weights)[~nans]
    #remove weights with zeros
    notzeros=np.nonzero(weights)
    values=values[notzeros]
    weights=weights[notzeros]
    #scale weights between 0 and 1
    weights=weights/np.sum(weights)
    if len(values)>0:
        return(weighted_median(values,weights))
    else:
        return(np.nan)


def insert_column_level(df,new_level,level_name):
    """
    inserts a new level to MultiIndex columns. Insertion of new index level possible by passing transposed dataframe and transposing again after.

    author: BZ
    
    Parameters
    ----------
    df : pd.DataFrame
        dataframe with multiindex columns.
    new_level : iterable 
        values for new level of length as df.columns
    level_name : str
        name for the new level.

    Returns
    -------
    pd.DataFrame

    """
    columns_df=pd.DataFrame(index=df.columns)
    columns_df[level_name]=new_level
    df.columns=pd.MultiIndex.from_frame(columns_df.reset_index())
    return(df)



def annotate_metadata_to_profiles(profiles, metadata_df, sample_col="Sample.ID"):
    """
    adds all metadata as multi-layer column names to profiles

    Parameters
    ----------
    profiles : pd.DataFrame
        DataFrame with IDs as rows, samples as columns and intensities as values.
    metadata_df : pd.DataFrame
        DataFrame with samples as rows and metadata variables as columns.
    sample_col : str , optional
        name of the column in metadata_df containing the sample names as in profiles.columns or profiles.columns level "sample". The default is "Sample.ID".

    Returns
    -------
    profiles : pd.DataFrame
        same format as input argument profiles with MultiIndex columns containing all variables from the metadata

    """
    if metadata_df.index.name is None:
        metadata_df.index.name="Index"
    metadata_index_orig=metadata_df.index.names
    metadata_df.reset_index(inplace=True)
    metadata_df.set_index(sample_col,inplace=True,drop=False)
    if "sample" in profiles.columns.names:
      profiles.columns=pd.MultiIndex.from_arrays(profiles.columns.to_frame().join(metadata_df.loc[profiles.columns.get_level_values("sample")],rsuffix="_metadata").T.values,names=profiles.columns.names + metadata_df.columns.tolist())
    else:
      profiles.columns=pd.MultiIndex.from_arrays(profiles.columns.to_frame().join(metadata_df.loc[profiles.columns],rsuffix="_metadata").T.values,names=profiles.columns.names + metadata_df.columns.tolist())
    metadata_df.set_index(metadata_index_orig,inplace=True)
    return(profiles)


def create_imputation_table(metadata_df,sample_hierarchies,missingness_thresholds="1/3;0.3;0.25_0.3",data_df=None,metadata_sample_column="Sample.Name",default_impu_method_ProtAbsence="knn",default_impu_method_MAR="knn",default_impu_method_MNAR="minval-noise",default_impu_method_ProtAbsenceMAR="knn",default_impu_method_ProtAbsenceMNAR="minval-noise",default_knn_k=15):
    """
    Create a template for an involved table to define different impuation methods across sample and precursor hierarchies.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Metadata dataframe mapping hierarchies to samples.
    sample_hierarchies : list
        list with names of columns in metadata_df to do hierarchical imputation across.
    missingness_thresholds : str
        str with thresholds for precursor missingness; ";" separates between sample_hierarchies, "_" separates within sample hierarchy. The default is "1/3;0.3;0.3_0.25".
    data_df : pd.DataFrame or None, optional
        wide format dataframe with IDs as rows and samples as columns, used to obtain samples in experiment; if None, all samples from metadata_df are taken into account. The default is None.
    metadata_sample_column : str, optional
        Name of the column containing sample names matching to columns in data_df. The default is "Sample.Name".
    default_impu_method_ProtAbsence : str, optional
        name of the imputation method used for missingness by 'ProtAbsence'. The default is 'knn'.
    default_impu_method_MAR : str, optional
        name of the imputation method used for missingness by 'MAR'. The default is 'knn'.
    default_impu_method_MNAR : str, optional
        name of the imputation method used for missingness by 'MNAR'. The default is 'minval-noise'.
    default_impu_method_ProtAbsenceMAR : str, optional
        name of the imputation method used for missingness by 'ProtAbsence/MAR'. The default is 'knn'.
    default_impu_method_ProtAbsenceMNAR : str, optional
        name of the imputation method used for missingness by 'ProtAbsence/MNAR'. The default is 'minval-noise'.
    default_knn_k : int, optional
        Default value for k used for kNN imputation.
        
    Returns
    -------
    pd.DataFrame.

    """
    
    if data_df is not None:
        samples = data_df.columns
    else:
        samples = metadata_df[metadata_sample_column]
    impu_table_df = pd.DataFrame(columns=sample_hierarchies + ["missingness_thresh","missingness_type","imputation_method","k"]) #TODO #add parameters for imputation methods
    #create dataframe storing classes (integers) for different missingness types
    missingness_types = ["ProtAbsence","MAR","MNAR","ProtAbsence/MAR","ProtAbsence/MNAR"]
    
    for s,sample_hierarchy in enumerate(sample_hierarchies):
        sample_hierarchy_values = metadata_df[metadata_df[metadata_sample_column].isin(samples)][sample_hierarchy].dropna().unique()
        sample_hierarchy_values = np.sort(sample_hierarchy_values)
        missingness_thresholds_s = [eval(x) for x in missingness_thresholds.split(";")[s].split("_")] # to convert e.g. 1/3 to float
        impu_table_part_df = pd.DataFrame(columns=impu_table_df.columns,index=range(len(sample_hierarchy_values)*len(missingness_types)*len(missingness_thresholds_s)))
        impu_table_part_df[sample_hierarchy] = list(itertools.chain.from_iterable([[x]*len(missingness_types)*len(missingness_thresholds_s) for x in sample_hierarchy_values]))
        impu_table_part_df["missingness_thresh"] = list(itertools.chain.from_iterable([[x]*len(missingness_types) for x in missingness_thresholds_s] * len(sample_hierarchy_values)))
        impu_table_part_df["missingness_type"] = missingness_types * len(missingness_thresholds_s) * len(sample_hierarchy_values)
        impu_table_df = pd.concat([impu_table_df,impu_table_part_df],axis=0,ignore_index=True)
    #set default parameter
    impu_table_df.loc[impu_table_df["missingness_type"]=="ProtAbsence","imputation_method"] = default_impu_method_ProtAbsence
    impu_table_df.loc[impu_table_df["missingness_type"]=="ProtAbsence/MAR","imputation_method"] = default_impu_method_ProtAbsenceMAR
    impu_table_df.loc[impu_table_df["missingness_type"]=="ProtAbsence/MNAR","imputation_method"] = default_impu_method_ProtAbsenceMNAR
    impu_table_df.loc[impu_table_df["missingness_type"]=="MAR","imputation_method"] = default_impu_method_MAR
    impu_table_df.loc[impu_table_df["missingness_type"]=="MNAR","imputation_method"] = default_impu_method_MNAR
    impu_table_df.loc[impu_table_df["imputation_method"]=="knn","k"] = default_knn_k
    return(impu_table_df)

def logger(log_file, module_output=str):
    """
    Iteratively callable function to add information to the existing logfile of the current preprocessing run
    :param log_file: str,
    path to existing logfile to append module information
    :param module_output: str,
    log_output from module to update the log file with module specific information like date, time, parameter values etc.
    """
    with open(log_file, "a") as logfile:  # open logfile
        logfile.write(module_output)  # append logfile
    return

def get_parser():
    """
    Generates the parser to parse arguments from the terminal to the Python Code

    Returns
    -------
    parser : object
        Parser containing input arguments for the initialization.
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--datapath", required=True, 
        help="Filepath to datafile that will be processed") 
    parser.add_argument(
        "-p", "--parampath", required=True, 
        help="Filepath to parameter json") 
    parser.add_argument(
        "-o", "--outpath", required=True, 
        help="Filepath for outputfile") 
    parser.add_argument(
        "-m", "--metapath", required=False, 
        help="Filepath to metadatafile") 
   # parser.add_argument(
   #     "-ma", "--maskpath", required=False,  #### CGenerate with datapath
   #     help="Filepath to maskfile") 
    parser.add_argument(
        "-l", "--log_file", required=True, 
        help="Filepath to logfile") 
    parser.add_argument(
        "-i", "--idmapping_file", required=False, 
        help="Filepath to ID mapping file") 
    parser.add_argument(
        "-c", "--counter", required=False, 
        help="Module counter", default=0)
    return parser

def get_parser_arguments(parser):
    """
    Extraction of parser arguments

    Parameters
    ----------
    parser : object
        Parser containing input arguments for the initialization.

    Returns
    -------
    data_df_fp : str
        File path to input data that should be processed.
    param_fp : str
        File path to parameter file
    out_dir_fp : str
        Directory for files generated during the module.
    meta_df_fp : str
        File path to the metadata file.
    log_file : str
        Log file.
    idmap_fp : str
        File path ID mapping file.
    counter : int
        Counter for the step, used to modfiy the seed and output file naming.
    outname : str
        Name of the output file: data_${counter}.tsv using a two digit counter.
    """
    
    args = parser.parse_args()
    data_df_fp = args.datapath.replace("\\","/")
    param_fp = args.parampath.replace("\\","/")
    out_dir_fp = args.outpath.replace("\\","/")
    try:
        meta_df_fp = args.metapath.replace("\\","/")
    except:
        meta_df_fp = ""
    try:
        idmap_fp = args.idmapping_file.replace("\\","/")
    except:
        idmap_fp = ""
    log_file = args.log_file.replace("\\","/")
    counter = args.counter
    outname = "data_"+"{:02d}".format(int(counter))+".tsv"
    return data_df_fp, param_fp, out_dir_fp, meta_df_fp, log_file, idmap_fp, int(counter), outname

def load_data(data_df_fp):
    """
    Loads data for processing.
    NOTE: At the moment only *.tsv is supported. This will be adapted in the future.

    Parameters
    ----------
    data_df_fp : str
        File path to input data that should be processed.

    Returns
    -------
    data_df : pd.DataFrame()
        Data Frame containing the data that are processed.
    """
    
    return pd.read_csv(data_df_fp,sep="\t",index_col=0)
    
def load_parameters(param_fp):
    """
    Read in parameter *.json file.

    Parameters
    ----------
    param_fp : str
        File path of the Parameter file.

    Returns
    -------
    parameters_dict : dict
        Dictionary containg parameters used for the model.
    """
    
    with open(param_fp) as f:
            parameters_dict = json.load(f)
    #### In python: 1==True and 0==False -> No need for conversion here!
    return parameters_dict
  
def load_id_mapping(idmap_fp):
    """
    Reads in ID mapping.

    Parameters
    ----------
    idmap_fp : str
        File path ID mapping file.

    Returns
    -------
    id_mapping_df : pd.DataFrame()
        Data frame for mapping precursors to proteins.

    """
    return pd.read_csv(idmap_fp,sep="\t")
    
def get_seed(parameters_dict, counter):
    """
    Generate the seed for the module.

    Parameters
    ----------
    parameters_dict : dict
        Dictionary containg parameters used for the model.
    counter : int
        Counter for the step, used to modfiy the seed and output file naming.

    Returns
    -------
    seed : int
        Seed used for the module calculated as parameters_dict["seed"] + counter.

    """
    return parameters_dict["seed"] + counter

def add_parameters_to_log(log_file, module, parameters_dict):
    """
    Writes parameters for the module to the log file

    Parameters
    ----------
    log_file : str
        Log file.
    module : str
        Module name.
    parameters_dict : dict
        Dictionary containg parameters used for the model.

    Returns
    -------
    module_output : str
        Resets module output after parameters are saved.

    """
    module_output = time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+module+" started.\n"
    module_output += "Parameters\n------------\n"
    for key,value in parameters_dict.items():
        module_output += str(key)+"\t\t\t"+str(value)+"\n"
    module_output += "------------\n"  
    
    logger(log_file, module_output=module_output)
    return ""

### Load Metadata without index
""""
def load_metadata(meta_df_fp, encoding="utf-8"):
    return pd.read_csv(meta_df_fp,sep="\t")
"""
### Load Mask
def load_mask(data_df_fp):
    """
    Loads the filtering/imputation mask.

    Parameters
    ----------
    data_df_fp : str
        File path to input data that should be processed. Name of the mask is generated as ${data_df_fp[:-4]}_mask.tsv.

    Returns
    -------
    masked_df : pd.DataFrame()
        DataFrame containing identifiers for masking steps.
    """
    
    return pd.read_csv(data_df_fp[:-4] + "_mask.tsv",sep="\t",index_col=0)

def load_mask_description(data_df_fp):
    """
    Loads description for the mask identifers

    Parameters
    ----------
    data_df_fp : str
        File path to input data that should be processed. Name of the mask description is generated as ${data_df_fp[:-4]}_mask_description.tsv.

    Returns
    -------
    mask_description_df : pd.DataFrame()
        Mapping of identifiers to corresponding operations.
    """
    
    return pd.read_csv(data_df_fp[:-4] + "_mask_description.tsv",sep="\t",index_col=0)
	
### Save mask
def save_mask(masked_df, out_dir_fp, outname):
    """
    Saves (updated) mask to disk,

    Parameters
    ----------
    masked_df : pd.DataFrame()
        DataFrame containing identifiers for masking steps.
    out_dir_fp : str
        Directory for files generated during the module.
    outname : str
        Name of the output file: ${outname}_mask.tsv using a two digit counter.

    Returns
    -------
    """
    
    masked_df.to_csv(os.path.join(out_dir_fp,outname[:-4] + "_mask.tsv"),sep="\t")
    return

def save_mask_description(mask_description_df, out_dir_fp, outname):
    """
    Saves (updated) mask description to disk,

    Parameters
    ----------
    mask_description_df : pd.DataFrame()
        Mapping of identifiers to corresponding operations.
    out_dir_fp : str
        Directory for files generated during the module.
    outname : str
        Name of the output file: ${outname}_mask_description.tsv using a two digit counter.

    Returns
    -------
    """
    
    mask_description_df.to_csv(os.path.join(out_dir_fp,outname[:-4] + "_mask_description.tsv"),sep="\t")
    return

def save_data(data_processed_df, out_dir_fp, outname):
    """
    Saves data after processing to disk. 
    Loads data for processing.
    NOTE: At the moment only *.tsv is supported. This will be adapted in the future.

    Parameters
    ----------
    data_processed_df : pd.DataFrame
        DataFrame containing the processed data.
    out_dir_fp : str
        Directory for files generated during the module.
    outname : str
        Name of the output file.

    Returns
    -------
    """
    
    data_processed_df.to_csv(os.path.join(out_dir_fp,outname),sep="\t")
    return

def init_module(module, get_metadata=False, get_id_mapping=False, get_mask=False, return_data_fp_only=False,
                skip_parser=False, data_df_fp=None, param_fp=None, out_dir_fp=None, meta_df_fp=None, log_file=None, idmap_fp=None, counter=None, outname=None):
    """
    Initialized the module and reads in all necessary data and arguments.

    Parameters
    ----------
    module : str
        Name of the module.
    get_metadata : bool, optional
        Switch to True to read-in metadata. The default is False.
    get_id_mapping : bool, optional
        Switch to True to read-in id_mapping. The default is False.
    get_mask : bool, optional
        Switch to True to read-in the mask. The default is False.
    return_data_fp_only : bool, optional
        Switch to True to read-in the mask. The default is False.
    skip_parser : bool, optional
        Switch to True to skip the parser and provide parser arguments manually
    data_df_fp : str, optional
        If skip_parser==True required. File path to input data that should be processed.
    param_fp : str, optional
        Directory containing parameter files. Parameter file is called as ${param_fp}/${out_dir_fp}_parameters.json
    out_If skip_parser==True required. dir_fp : str, optional
        If skip_parser==True required. Directory for files generated during the module.
    meta_df_fp : str, optional
        If skip_parser==True and load_metadata==True required. File path to the metadata file.
    log_file : str, optional
        If skip_parser==True required. Log file.
    idmap_fp : str, optional
        If skip_parser==True and load_id_mapping ==True required. File path ID mapping file.
    counter : int, optional
        If skip_parser==True required. Counter for the step, used to modfiy the seed and output file naming.
    outname : str, optional
        If skip_parser==True required. Name of the output file: data_${counter}.tsv using a two digit counter.

    Returns
    -------
    data_df : pd.DataFrame()
        Data Frame containing the data that are processed. Only returned if return_data_fp==False.
    parameters_dict : dict
        Dictionary containg parameters used for the model.
    id_mapping_df : pd.DataFrame()
        Data frame for mapping precursors to proteins.
    seed : int
        Seed used for the module calculated as parameters_dict["seed"] + counter.
    module_output : str
        Resets module output after parameters are saved.
    metadata_df : pd.DataFrame()
        DataFrame containing the metadata.
    masked_df : pd.DataFrame()
        DataFrame containing identifiers for masking steps.
    mask_description_df : pd.DataFrame()
        Mapping of identifiers to corresponding operations.
    out_dir_fp : str
        Directory for files generated during the module.
    outname : str
        Name of the output file: data_${counter}.tsv using a two digit counter.
    log_file : str
        Log file.
    counter : int
        Counter for the step, used to modfiy the seed and output file naming.
    data_df_fp : str
        File path to input data that should be processed. Only returned in module_filtering_long.py (return_data_fp==True).
    """
    
    # Parser
    if skip_parser:
        if not data_df_fp:
            raise ValueError("###ERROR###: data_df_fp required.")
        if not param_fp:
            raise ValueError("###ERROR###: param_fp required.")
        if not out_dir_fp:
            raise ValueError("###ERROR###: out_dir_fp required.")
        if not meta_df_fp and get_metadata:
            raise ValueError("###ERROR###: meta_df_fp required.")
        if not log_file:
            raise ValueError("###ERROR###: log_file required.")
        if not idmap_fp and get_id_mapping:
            raise ValueError("###ERROR###: idmap_fp required.")
        if not counter:
            raise ValueError("###ERROR###: counter required.")
        if not outname:
            raise ValueError("###ERROR###: outname required.")
    else:
       parser = get_parser() 
       data_df_fp, param_fp, out_dir_fp, meta_df_fp, log_file, idmap_fp, counter, outname = get_parser_arguments(parser)
    # Load data
    if not return_data_fp_only:
        data_df = load_data(data_df_fp)
    else:
        data_df = None
    # Load parameters
    parameters_dict = load_parameters(param_fp)
    # Load ID mapping
    if get_id_mapping:
        id_mapping_df = load_id_mapping(idmap_fp)
    else:
        id_mapping_df = []
    # Get seed
    seed = get_seed(parameters_dict, counter)
    # Write Parameters to log-file
    module_output = add_parameters_to_log(log_file, module, parameters_dict)
    # Load meta data
    if get_metadata:
        metadata_df, log = load_metadata(metadata_fp=meta_df_fp,
                                         metadata_sheet=parameters_dict["metadata_sheet"], 
                                         filename_col=parameters_dict["metadata_filename_column"], 
                                         encoding=parameters_dict["metadata_encoding"])
        module_output += log
    else:
        metadata_df = []
    # Load mask
    if get_mask:
        masked_df = load_mask(data_df_fp)
        mask_description_df = load_mask_description(data_df_fp)
    else:
        masked_df = []
        mask_description_df = []
    if return_data_fp_only:
        return data_df_fp, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, mask_description_df, out_dir_fp, outname, log_file, counter
    else:
        return data_df, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, mask_description_df, out_dir_fp, outname, log_file, counter
	    
	    
def close_module(data_processed_df, masked_df, mask_description_df, log_file, module_output, out_dir_fp, outname, mask=False, exit_module=True):
    """
    Saves data to disk and exits the module.

    Parameters
    ----------
    data_processed_df : pd.DataFrame()
        Data Frame containing processed data that are processed.
    masked_df : pd.DataFrame()
        DataFrame containing identifiers for masking steps.
    mask_description_df : pd.DataFrame()
        Mapping of identifiers to corresponding operations.
    log_file : str
        Log file.
    module_output : str
        Module output to be written to the log file.
    out_dir_fp : str
        Directory for files generated during the module.
    outname : str
        Name of the output file: data_${counter}.tsv using a two digit counter.
    mask : bool, optional
        Save the mask to disk. The default is False.
	exit_module : bool
		call sys.exit, switch to False if run within Python. The default is True.

    Returns
    -------
    """
    
    # Update log file
    logger(log_file, module_output=module_output)
    # Save data
    save_data(data_processed_df, out_dir_fp, outname)
    # Save mask
    if mask:
        save_mask(masked_df, out_dir_fp, outname)
        save_mask_description(mask_description_df, out_dir_fp, outname)
    # Close function
    logger(log_file, time.strftime("\n%Y-%m-%d%t%H:%M:%S") + "\tModule finished as expected.\n")
    if exit_module:
        sys.exit(0)
