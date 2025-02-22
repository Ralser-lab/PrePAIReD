"""
Input: X
Output: X

functions:
    - create_analysis_folder: if not existent, creates a folder for analysis
    - replace_invalid_characters: replaces characters that cannot be decoded in standard encoding
    - create_logfile: creates a textfile documenting the filtering and pre-processing analysis
    - create_filtering_matrix: creates a wide-format .tsv file which keeps track for every ID in every sample if it has been measured, filtered (and why) or imputed
    - rename_filter_mask_columns: renames filenames in filtering_matrix to shorter sample names as in the filtered output
    - create_filename_annotations_file: if no metadata present, creates a metadata .tsv file containing the original filenames, shorter sample names and sample groups
    - check_filename_annotations_on_unique_indices: sample and group names in created metadata can be assigned manually (e.g. via Excel). This function checks whether sample names are unique
    - load_metadata: loads existing metadata, either automatically generated tsv-file or manually curated xlsx-file

"""

import re
import numpy as np
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import pandas as pd
import numpy as np
import os.path
import os
import glob


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
    df = base**df
    return(df)

def mad(x):
    """
    Xalculate median absolute deviation of a vector x
    
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
def create_filtering_matrix(file,metadata_df,min_peptides,values,ids=["Precursor.Id","Stripped.Sequence","Protein.Group","Genes"],filename_col="filename",diann_sample_column="File.Name",logfilename=None):
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
        minimum number of precursors required to quantify a protein.
    values : str
        values extracted from DIA-NN report.tsv, e.g. Precursor.Quantity.
    ids : list, optional
        list of ids extracted from DIA-NN report.tsv. The default is ["Precursor.Id","Stripped.Sequence","Protein.Group","Genes"].
    filename_col : str, optional
        column name with complete filenames in metadata. The default is "filename".
    diann_sample_column : str, optional
        name of the sample names in the DIA-NN run to match with metadata, either "Run" or "File.Name". The default is "File.Name".
    logfilename : str
        name (path + filename) of the logfile. The default is None, disabling logging.
        
    Returns
    -------
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters

    """
    print("Creating filtering matrix...")
    #get sorted File.Names
    samples=[x.replace("\\","/").split("/")[-1] for x in file[diann_sample_column].unique()]
    if len(set(samples)) > len(set(metadata_df.index)):
        print("WARNING: Not all samples in DIA-NN output annotated in metadata")
    elif len(set(samples)) < len(set(metadata_df.index)):
        print("WARNING: Not all samples in metadata found in DIA-NN output")
    #get intersection between filenames and metadata
    samples=list(set(samples) & set(metadata_df.index))
    #sort File.Names
    samples.sort()
    columns_filenames=metadata_df.loc[samples,filename_col]
    mask = pd.DataFrame(index=file.set_index(ids).index.unique(), columns=columns_filenames,data=0)
    # dataframe explaining reasons why peptide was excluded
    mask_description = pd.DataFrame(columns=["filter","value"])
    #mask_description["filter"] = ["not measured","measured and kept", "excessing_NaNs","Q.Value","Global.Q.Value","GG.Q.Value","PG.Q.Value","Protein.Q.Value","Global.PG.Q.Value","Lib.Q.Value","Lib.PG.Q.Value","Precursor.Charge","Proteotypic","Min. "+str(min_peptides)+" Peptides","Zero "+values,"Imputed"]
    #reordered to same order as filtering steps
    mask_description["filter"] = ["not measured","measured and kept","Proteotypic","Q.Value","Global.Q.Value","GG.Q.Value","PG.Q.Value","Protein.Q.Value","Global.PG.Q.Value","Lib.Q.Value","Lib.PG.Q.Value","PEP","Ms1.Profile.Corr","Mass.Evidence","Ms1.Profile.Corr + Mass.Evidence","Zero "+values,"Precursor.Charge","excessing_NaNs","Min. "+str(min_peptides)+" Peptides","missing not at random","Imputed"]
    mask_description["value"]=range(0,len(mask_description))
    #write mask description into logfile
    if logfilename is not None:
        logfile=open(logfilename,"a")
        logfile.write("\nFiltering matrix to keep track over excluded precursors and proteins created.")
        logfile.write("\nReasons for exclusion of precursors as specified in filter mask\n")
        logfile.write(mask_description.to_string(index=False))
        logfile.write("\n-----------------------------------------------\n\n")
        logfile.close()
    mask_description.set_index("filter",inplace=True)
    mask_description.loc["Imputed","value"]=99
    return([mask,mask_description])


def rename_filter_mask_columns(filter_mask,metadata_df,filename_col="filename",sample_col="sample",sample_group_col="sample group"):
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
    sample_group_col: str, optional
        column name with sample group names in metadata. The default is "sample group". 

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
    filter_mask.columns=[metadata_df.set_index(filename_col,drop=False).loc[filter_mask.columns,sample_col],metadata_df.set_index(filename_col,drop=False).loc[filter_mask.columns,sample_group_col]]
    filter_mask.columns.names=["sample","group"] ### added 
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


#load metadata file and assign columns as present in metadata_df (output from create_filename_annotations_file)
# metafilename is name including path of metadata file
# if metadata file is an excel file, metadatasheet is the excel sheet in the file providing the annotations
# ending is file-ending of files in dia-nn output, e.g. ".d", ".rawIdx.wiff.dia"
def load_metadata(metafilename=None,metadata_sheet="Sheet1",ending="",sample_column="Sample Name",group_column="PatID",label_column=None, filename_col="filename", encoding="utf-8", logfilename=None):
    """
    load the metadata file
    
    author: BZ
    
    Parameters
    ----------
    metafilename : str, optional
        name of the metadata file. The default is None, if None a GUI opens to select the file.
    metadata_sheet : str, optional
        if metadata is an excel file, metadata_sheet specifies the sheetname. The default is "Sheet1".
    ending : str, optional
        file names in metadata may have capped file endings compared to DIA-NN output report, ending specifies the file ending to match report.tsv["File.Name"]. The default is "".
    group_column : str, optional
        name of column in metadata specifying the sample group, used for e.g. separate plots per group. The default is "PatID".
    label_column : str, optional
        name of column in metadata specifying a label, used for e.g. coloring in plots. The default is None.
    filename_col : str, optional
        column name with complete filenames in metadata. The default is "filename".
    encoding : str, optional
        encoding passed to pandas read_csv
    logfilename : str
        name (path + filename) of the logfile. The default is None, disabling logging.
        
    Returns
    -------
    metafilename : str, optional
        name of the metadata file.
    metadata_df : pd.DataFrame
        the metadata file as dataframe

    """
    if metafilename is None:
        print("Please select file containing meta data.")
        root = Tk()
        root.withdraw()
        metafilename = askopenfilename(
            initialdir="S:/AG/AG-CF-HTMS/CoreFacility/")
        root.destroy()
    # read tsv file
    if metafilename[-4:] in [".txt",".tsv"]:
        metadata_df = pd.read_csv(metafilename, sep="\t", encoding=encoding)
    elif metafilename[-4:]==".csv":
        metadata_df = pd.read_csv(metafilename, sep=",", encoding=encoding)
    #read xlsx file
    elif metafilename[-5:]==".xlsx":
        metadata_excel=pd.ExcelFile(metafilename)
        metadata_df=metadata_excel.parse(metadata_sheet)
        metadata_excel.close()
        #add columns
        if "Path" in metadata_df.columns and "File.Name" in metadata_df.columns:
            metadata_df[filename_col]=metadata_df["Path"]+metadata_df["File.Name"]+ending
        '''
        #add sample name, e.g. ??
        if sample_column is not None:
            metadata_df["name on x-axis"]=metadata_df[sample_column]
        #add group for separate plots in QC
        if group_column is not None:
            metadata_df["sample group"]=metadata_df[group_column]
        '''
        #add label, for colouring
        if label_column is not None:
            metadata_df["label"]=metadata_df[label_column]
    else:
        print("Please implement how how to read metadata of filetype "+metafilename[-5:])
    #exclude all annotations without filename (e.g. not measured yet)
    if filename_col in metadata_df.columns:
        metadata_df=metadata_df.dropna(subset=[filename_col])
    if "FileName_trimmed" not in metadata_df.columns and filename_col in metadata_df.columns:
        metadata_df["FileName_trimmed"] = list(map(lambda x: x.replace("/", "\\").split("\\")[-1], metadata_df[filename_col].tolist()))
    if "FileName_trimmed" in metadata_df.columns:
        metadata_df.set_index("FileName_trimmed", inplace=True)
    # drop duplicated indices
    if sum(metadata_df.index.duplicated())>0:
        print(str(sum(metadata_df.index.duplicated()))+" indices duplicated in metadata. Duplicated indices are being dropped.")
        if logfilename is not None:
            logfile=open(logfilename,"a")
            logfile.write(str(sum(metadata_df.index.duplicated()))+" indices duplicated in metadata. Duplicated indices are being dropped.\n")
            logfile.close()
        metadata_df=metadata_df.drop_duplicates() 
    print("Meta-data/Profile Annotations loaded successfully.")
    return([metafilename,metadata_df])

