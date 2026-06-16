#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:22:54 2026

@author: oliverlemke
"""

import pandas as pd
import pickle
import os
import numpy as np
from pyarrow.parquet import ParquetFile
import pyarrow as pa

def load_metadata(metadata_fp, metadata_sheet="Sheet1", encoding="utf-8"):
    """
    Loads metadata.
    
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

    """
    # read tsv file
    try:
        if metadata_fp[-4:] in [".txt",".tsv"]:
            metadata_df = pd.read_csv(metadata_fp, sep="\t", encoding=encoding)
        elif metadata_fp[-4:]==".csv":
            metadata_df = pd.read_csv(metadata_fp, sep=",", encoding=encoding)
        #read xlsx file
        elif metadata_fp[-5:]==".xlsx":
            metadata_excel=pd.ExcelFile(metadata_fp)
            #get sheets, check if sheet in there
            if metadata_sheet not in metadata_excel.sheet_names:
                return None
            metadata_df=metadata_excel.parse(metadata_sheet)
            metadata_excel.close()
        else:
            return None
    except:
        return None
    #exclude all annotations without filename (e.g. not measured yet)
    # drop duplicated indices
    if sum(metadata_df.index.duplicated())>0:
        metadata_df=metadata_df.drop_duplicates()
    return metadata_df

def load_data_header(ms_data_fp):
    """
    Loads data header of MS-file.

    Parameters
    ----------
    ms_data_fp : str
        File path to input data that should be processed.

    Returns
    -------
    ms_data_columns : pd.DataFrame()
        Data Frame containing the data that are processed.
        
    """
    
    if ms_data_fp.endswith(".tsv"):
        columns = pd.read_csv(ms_data_fp, sep="\t", nrows=0).columns.tolist()
    elif ms_data_fp.endswith(".parquet"):  # read parquet files
        pf = ParquetFile(ms_data_fp)
        first_row = next(pf.iter_batches(batch_size=1))
        df = pa.Table.from_batches([first_row]).to_pandas()
        columns = df.columns.tolist()
    return columns
    
def load_parameter_dictionaries(param_path_fp):
    """
    Loads parameter dictionaries.
    
    Parameters
    ----------
    param_path_fp : str
        File path to parameter dictionaries.

    Returns
    -------
    flag_dict : dict
        Flag dictionary.
    parameters_dict : dict
        parameters dictionary.
    
    """
    
    flag_dict = pickle.load(open(os.path.join(param_path_fp,"dict_flags.pkl"),"rb"))
    flag_dict["FILTERING_WIDE"]["filtering_wide_methods"]["Default"] = None
    parameters_dict = pickle.load(open(os.path.join(param_path_fp,"dict_parameters.pkl"),"rb"))
    return flag_dict, parameters_dict

def get_steps(workflow_data, params_dict, rerun_counter=0):
    """
    Extract the steps of the workflow. Summarises certain steps for less complexity.
    
    Parameters
    ----------
    workflow_data : list
        Workflow setup.
    params_dict: dict
        Parameter dictionary.
    rerun_counter : int
        Last counter of the prior run. The default is 0.
        
    Returns
    -------
    steps : dict
        Dictionary containing all steps with an identifier.
    precursor_metadata: dict
        Metadata for filtering precursors by metadata. 
    language : dict
        Programming language of the selected method.
    methods: dict
        Collection of methods.
    step_collection : dict
        Collection of counters.
        
    """
    steps = {}
    fw_methods = []
    #collection = []
    precursor_metadata_col = ""
    precursor_metadata_opr = ""
    precursor_metadata_val = ""
    precursor_metadata = {}
    language = {}
    methods = {}
    step_collection = {}
    
    counter_step = 0 + rerun_counter
    
    for ind,item in enumerate(workflow_data):
        step_collection[ind] = counter_step
        step = item["name"]
        if step != "Filtering_wide":
            if fw_methods:
                steps.update({counter_step:"Filtering_wide"})
                precursor_metadata.update({counter_step:{"col":precursor_metadata_col.strip(";"),"opr":precursor_metadata_opr.strip(";"),"val":precursor_metadata_val.strip(";")}})
                language.update({counter_step:"Python"})
                methods.update({counter_step:fw_methods})
                
                fw_methods = []
                precursor_metadata_col = ""
                precursor_metadata_opr = ""
                precursor_metadata_val = ""
                counter_step +=1
            
            step_collection[ind] = counter_step
            steps.update({counter_step:step})
            counter_step +=1
            
        elif item["params"]["filtering_wide_methods"] not in fw_methods or item["params"]["filtering_wide_methods"] == "filter_precursors_by_metadata":
            fw_methods.append(item["params"]["filtering_wide_methods"])
            if item["params"]["filtering_wide_methods"] == "filter_precursors_by_metadata":
                family = item["params"]["precursor_filter_col"]
                if family[0]["parent"] and family[0]["child"]and family[0]["child"]:
                    precursor_metadata_col += family[0]["parent"]+";"
                    precursor_metadata_opr += family[0]["child"]+";"
                    precursor_metadata_val += family[1]["child"]+";"
        else:
            steps.update({counter_step:step})
            precursor_metadata.update({counter_step:{"col":precursor_metadata_col.strip(";"),"opr":precursor_metadata_opr.strip(";"),"val":precursor_metadata_val.strip(";")}})
            methods.update({counter_step:fw_methods})
            fw_methods = []
            precursor_metadata_col = ""
            precursor_metadata_opr = ""
            precursor_metadata_val = ""
            
            fw_methods.append(item["params"]["filtering_wide_methods"])
            if item["params"]["filtering_wide_methods"] == "filter_precursors_by_metadata":
                family = item["params"]["precursor_filter_col"]
                if family[0]["parent"] and family[0]["child"]and family[0]["child"]:
                    precursor_metadata_col += family[0]["parent"]+";"
                    precursor_metadata_opr += family[0]["child"]+";"
                    precursor_metadata_val += family[1]["child"]+";"
            
            counter_step +=1
            step_collection[ind] = counter_step
            
        if step == "Filtering_long":
            language.update({counter_step-1:"Python"})
        elif step == "Filtering_wide":
            language.update({counter_step-1:"Python"})
        elif step == "Normalisation":
            language.update({counter_step-1:params_dict["normalization_method"][item["params"]["normalization_method"]]["Language"]})
            methods.update({counter_step-1:item["params"]["normalization_method"]})
        elif step == "Imputation":
            language.update({counter_step-1:params_dict["imputation_method"][item["params"]["imputation_method"]]["Language"]})
            if item["params"]["simple_imputation"]=="0":
                methods.update({counter_step-1:[item["params"]["imputation_method"],item["params"]["eval_missingness_by_detection_threshold"]]})
            else:
                methods.update({counter_step-1:[item["params"]["imputation_method"]]})
        elif step == "Correction":
            language.update({counter_step-1:params_dict["correct_method"][item["params"]["correct_method"]]["Language"]})
            methods.update({counter_step-1:item["params"]["correct_method"]})
        elif step == "Summarisation":
            language.update({counter_step-1:params_dict["summarization_method"][item["params"]["summarization_method"]]["Language"]})
            methods.update({counter_step-1:item["params"]["summarization_method"]})
    
    if fw_methods:
        steps.update({counter_step:"Filtering_wide"})
        language.update({counter_step:"Python"})
        methods.update({counter_step:fw_methods})
    
    if precursor_metadata_col:
        precursor_metadata.update({counter_step:{"col":precursor_metadata_col.strip(";"),"opr":precursor_metadata_opr.strip(";"),"val":precursor_metadata_val.strip(";")}})
    
    return steps, precursor_metadata, language, methods, step_collection

def write_config(workdir, steps, language, rerun=False, rerun_start=""):
    """
    Writes configuration file to disk.
    
    Parameters
    ----------
    workdir : str
        Working directory.
    steps: dict
        Collection of steps.
    language: dict
        Programming language for each step.
    rerun : bool
        If True, the configuration is written for a rerun. The default is False
    rerun_start : str
        Step of prior run from which rerun should be started.
        
    Returns
    -------
    step_names : dict
        Dictionary containg the full step names.
        
    """
    
    step_names = {}
    with open(os.path.join(workdir,"Config.txt"),"w") as file:
        file.write("Run = ")
        if rerun:
            file.write(rerun_start+"\n")
        else:
            file.write("new\n")
        file.write("order = ")
        for ind,(key,step) in enumerate(steps.items()):
            if step != "General":
                step_name = "{:02d}".format(key)
                step_name += "_"
                step_name += step.upper()
                
                step_names.update({key:step_name})
                
                file.write(step_name)
                file.write("_")
                if language[key] == "Python":
                    file.write("py")
                elif language[key] == "R":
                    file.write("R")
                if ind<len(steps)-1:
                    file.write(", ")
    return step_names

def write_evaluation_config(workdir, step_names, methods, rerun=False, rerun_start=None, docker_setup=False):
    """
    Writes evaluation configuration to disk.
    
    Parameters
    ----------
    workdir : str
        Working directory.
    step_names : dict
        Dictionary containg the full step names.
    methods : dict
        Dictionary containg all methods.
    erun : bool
        If True, the configuration is written for a rerun. The default is False.
    rerun_start : str
        Step of prior run from which rerun should be started.
    docker_setup : bool
        If True, pathes are adjusted to match NextFlow setup. The default is False.
    
    Returns
    -------
    None
    
    """
    
    if docker_setup:
        outpath = ""
    else:
        outpath = workdir
    
    with open(os.path.join(workdir,"Config_evaluation.txt"),"w") as file:
        file.write("%%%%%%%%%%00_GENERAL%%%%%%%%%%"+"\n")
        file.write("%-%-%-%"+"\n")
        file.write("param_fp = "+os.path.join(outpath,"Parameter_Files","00_GENERAL_Parameters.json"+"\n"))
    
        steps_ = []
        for ind,(key,step) in enumerate(step_names.items()):
            if docker_setup:
                steps_.append("")
            else:
                steps_.append(step)
            file.write("%%%%%%%%%%"+step+"%%%%%%%%%%"+"\n")
            file.write("%-%-%-%"+"\n")
            if ind == 0:
                if rerun:
                    if docker_setup:
                        file.write("data_origin = "+os.path.join(outpath,"","data_"+rerun_start.split("_")[0]+".tsv")+"\n")
                    else:
                        file.write("data_origin = "+os.path.join(outpath,rerun_start,"data_"+rerun_start.split("_")[0]+".tsv")+"\n")
                else:
                    file.write("data_origin = None"+"\n")
            else:
                file.write("data_origin = "+os.path.join(outpath,steps_[ind-1],"data_{:02d}".format(key-1)+".tsv")+"\n")
            file.write("data_processed = "+os.path.join(outpath,steps_[ind],"data_{:02d}".format(key)+".tsv")+"\n")
            file.write("param_fp = "+os.path.join(outpath,"Parameter_Files",step+"_Parameters.json")+"\n")
            file.write("%-%-%-%"+"\n")
            
            module = step.split("_",1)[1]
            
            # F01
            if module in ["FILTERING_LONG","FILTERING_WIDE","NORMALISATION","IMPUTATION","CORRECTION","SUMMARISATION"]:
                file.write("plot_int_box = 1"+"\n")
                file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F01.png")+"\n")
                if module == "FILTERING_LONG":
                    file.write("size = (5,10)"+"\n")
                else:
                    file.write("size = (10,10)"+"\n")
            else:
                file.write("plot_int_box = 0"+"\n")
                file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F01.png")+"\n")
            
            file.write("%-%-%-%"+"\n")
            
            # F02
            if module in ["FILTERING_LONG","FILTERING_WIDE"]:
                file.write("plot_filter_mask = 1"+"\n")
            else:
                file.write("plot_filter_mask = 0"+"\n")
            file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F02.png")+"\n")
            
            file.write("%-%-%-%"+"\n")
            
            # F03
            if module in ["NORMALISATION","IMPUTATION","CORRECTION","SUMMARISATION"]:
                file.write("plot_cv = 1"+"\n")
            else:
                file.write("plot_cv = 0"+"\n")
            file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F03.png"+"\n"))
    
            file.write("%-%-%-%"+"\n")
            
            # F04
            if module in ["NORMALISATION","IMPUTATION","CORRECTION","SUMMARISATION"]:
                file.write("plot_pca = 1"+"\n")
                file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F04.png"+"\n"))
                if module == "NORMALISATION":
                    file.write("continuous = True"+"\n")
                    file.write("coloring_property = Missingness"+"\n")
                    file.write("use_metadata_coloring = False"+"\n")
                elif module == "IMPUTATION":
                    file.write("continuous = True"+"\n")
                    file.write("coloring_property = Mean_Intensities"+"\n")
                    file.write("use_metadata_coloring = False"+"\n")
                elif module == "CORRECTION":
                    file.write("continuous = None"+"\n")
                    file.write("coloring_property = None"+"\n")
                    file.write("use_metadata_coloring = True"+"\n")
                elif module == "SUMMARISATION":
                    file.write("continuous = False"+"\n")
                    file.write("coloring_property = None"+"\n")
                    file.write("use_metadata_coloring = True"+"\n")
            else:
                file.write("plot_pca = 0"+"\n")
                file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F04.png"+"\n"))
                file.write("continuous = True"+"\n")
                file.write("coloring_property = Missingness"+"\n")
                file.write("use_metadata_coloring = False"+"\n")
                
            file.write("%-%-%-%"+"\n")
            
            # F05
            if module in ["NORMALISATION"]:
                file.write("plot_distribution = 1"+"\n")
            else:
                file.write("plot_distribution = 0"+"\n")
            file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F05.png"+"\n"))
            
            file.write("%-%-%-%"+"\n")
            
            # F06
            if module in ["NORMALISATION","CORRECTION"]:
                if module == "NORMALISATION":
                    file.write("plot_scatter_series = 1"+"\n")
                else:
                    if methods[key] in ["adjust_confounders","adjust_confounders_batchwise"]:
                        file.write("plot_scatter_series = 1"+"\n")
                    else:
                        file.write("plot_scatter_series = 0"+"\n")
            else:
                file.write("plot_scatter_series = 0"+"\n")     
            file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F06.png"+"\n"))
            file.write("plot_by = Correlation"+"\n")
            
            file.write("%-%-%-%"+"\n")
            
            # F07
            if module in ["IMPUTATION"]: 
                if "1" in methods[key]:
                    file.write("plot_decision_boundary = 1"+"\n")
                else:
                    file.write("plot_decision_boundary = 0"+"\n")
            else:
                file.write("plot_decision_boundary = 0"+"\n")
            file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F07.png"+"\n"))
            
            file.write("%-%-%-%"+"\n")
            
            # F08
            if module in ["CORRECTION"]:
                if methods[key] in ["adjust_confounders","adjust_confounders_batchwise"]:
                    file.write("plot_correlation = 1"+"\n")
                else:
                    file.write("plot_correlation = 0"+"\n")
            else:
                file.write("plot_correlation = 0"+"\n")
            file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F08.png"+"\n"))
            
            file.write("%-%-%-%"+"\n")
            
            # F09
            if module in ["NORMALISATION"]: 
                if methods[key] == "median_scaling":
                    file.write("plot_distribution_hist = 1"+"\n")
                    file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F09.png"+"\n"))
                    file.write("normalization_factors = "+os.path.join(outpath,steps_[ind],"Normalization_factors.tsv"+"\n"))
                else:
                    file.write("plot_distribution_hist = 0"+"\n")
                    file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F09.png"+"\n"))
            else:
                file.write("plot_distribution_hist = 0"+"\n")
                file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F09.png"+"\n"))
                
            file.write("%-%-%-%"+"\n")
            
            # F10
            if module in ["FILTERING_LONG","FILTERING_WIDE"]: 
                file.write("plot_abundance_vs_completeness = 1"+"\n")
                file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F10.png"+"\n"))
                if module == "FILTERING_LONG":
                    file.write("gridsize = 50"+"\n")
                    file.write("n_bins = 50"+"\n")
                elif module == "FILTERING_WIDE":
                    file.write("gridsize = 25"+"\n")
                    file.write("n_bins = 25"+"\n")
            else:
                file.write("plot_abundance_vs_completeness = 0"+"\n") 
                file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F10.png"+"\n"))
                
            file.write("%-%-%-%"+"\n")
            
            # F11
            if module in ["FILTERING_LONG","FILTERING_WIDE"]: 
                file.write("plot_overall_completeness = 1"+"\n")
            else:
                file.write("plot_overall_completeness = 0"+"\n") 
            file.write("outfile = "+ os.path.join(outpath,steps_[ind],"report_folder","F11.png"+"\n"))            
    return

def write_parameter_file(workdir, workflow_data, step_names, step_collection, methods, precursor_metadata, flags_dict, params_dict):
    """
    Writes parameter file to disk.
    
    Parameters
    ----------
    workdir : str
        Working directory.
    workflow_data : list
        Workflow setup.
    step_names : dict
        Dictionary containg the full step names.
    step_collection : dict
        Collection of counters.
    methods : dict
        Dictionary containg all methods.
    precursor_metadata: dict
        Metadata for filtering precursors by metadata.  
    flags_dict : dict
        Dictionary for flag variables.
    params_dict : dict
        Dictionary for parameters.
    
    Returns
    -------
    None
    
    """
    
    simple_imputation_dropouts = ["imputation_method_HA","imputation_method_LA","imputation_method_ProtAbsence","imputation_method_ProtAbsenceHA","imputation_method_ProtAbsenceLA","eval_missingness_on_protein_presence","eval_missingness_by_detection_threshold","hierarchical_imputation","remember_missingness_type_in_mask","level_for_intensity_dependent_missingness_estimation","save_intermediate_output_for_QA","impute_over_whole_dataset","missingness_type_protein_absence","define_protein_missingness_on_filter_mask","intensity_dependent_missingness_decision_boundary","intensity_dependent_missingness_detection_proportion_cap","intensity_dependent_missingness_fit_method","intensity_dependent_missingness_n_knots"]
    step_ = None
    skip_fw_method = False
    skip_prec_metadata = False
    
    with open(os.path.join(workdir,"Parameters.txt"),"w") as file:
        for key,item in enumerate(workflow_data):
            # Reset checks
            skip_simple = False
            step = item["name"].upper()
            
            if step == "GENERAL":
                file.write("%%%%%%%%%%00_GENERAL%%%%%%%%%%"+"\n")
            elif step_ != step_names[step_collection[key]]:
                file.write("\n")
                file.write("%%%%%%%%%%"+step_names[step_collection[key]]+"%%%%%%%%%%"+"\n")
                skip_fw_method = False
                skip_prec_metadata = False
            if step != "GENERAL":
                step_ = step_names[step_collection[key]]    
            
            for param, el in item["params"].items():
                
                if param in flags_dict[step]:
                    current_flag = param
                    current_flag_value = el
                    flag = True
                    process = True
                    if current_flag == "simple_imputation" and current_flag_value == "1":
                        skip_simple = True
                    if skip_simple and param in simple_imputation_dropouts:
                        process = False
                    if current_flag == "filtering_wide_methods" and not skip_fw_method:
                        file.write("filtering_wide_methods = ")
                        for ind,method in enumerate(methods[step_collection[key]]):
                            file.write(method)
                            if ind < len(methods[step_collection[key]])-1:
                                file.write(";")
                        file.write("\n")
                        skip_fw_method = True
                    if skip_fw_method:
                        process = False
                        
                else:
                    flag = False
                    if skip_simple and param in simple_imputation_dropouts:
                        process = False
                    elif param == "precursor_filter_col" and not skip_prec_metadata and "filter_precursors_by_metadata" in methods[step_collection[key]]:
                        file.write("precursor_filter_col = "+precursor_metadata[step_collection[key]]["col"]+"\n")
                        file.write("filter_operator = "+precursor_metadata[step_collection[key]]["opr"]+"\n")
                        file.write("precursor_filter_value = "+precursor_metadata[step_collection[key]]["val"]+"\n")
                        skip_prec_metadata = True
                        process = False
                    elif param in ["precursor_filter_col", "filter_operator", "precursor_filter_value"] and skip_prec_metadata:
                        process = False
                    elif param in params_dict[current_flag][current_flag_value]:
                        process = True
                    else:
                        process = False
                
                if process:
                    file.write(param)
                    file.write(" = ")
                    # List and dictionary handling
                    if type(el) != list and type(el) != tuple:
                        file.write(str(el))
                        file.write("\n")
                    else:
                        # Dictionary handling
                        if type(el[0]) == dict:
                            if "children" in el[0]:
                                parent_vals = ""
                                child_vals = {i:"" for i in range(len(el[0]["children"]))}
                                for el_ in el:
                                    parent_vals += str(el_["parent"])+";"
                                    for i in range(len(el_["children"])):
                                        child_vals[i] += str(el_["children"][i]["child"])+";"
                                
                                el = [{"parent":parent_vals.strip(";"), "child":child_vals[i].strip(";")} for i in child_vals]
                            if el[0]["parent"] and el[0]["child"]:
                                if flag:
                                    var_type = flags_dict[step][param]["Type"]
                                else:
                                    var_type = params_dict[current_flag][current_flag_value][param]["Type"]
                                if "str" in var_type:
                                    file.write(el[0]["parent"])
                                elif "None" in var_type and (el[0]["parent"] is None or el[0]["parent"] == "None"):
                                    file.write(el[0]["parent"])
                                else:
                                    if ";" in el[0]["parent"]:
                                        file.write(str(el[0]["parent"].split(";")))
                                    else:
                                        file.write(str([el[0]["parent"]]))
                            else:
                                file.write("None")
                            file.write("\n")
                            if flag:
                                children = flags_dict[step][param]["Child"]
                            else:
                                children = params_dict[current_flag][current_flag_value][param]["Child"]
                            for child,val in zip(children,el):
                                if flag:
                                    var_type = flags_dict[step][child]["Type"]
                                else:
                                    var_type = params_dict[current_flag][current_flag_value][child]["Type"]
                                file.write(child)
                                file.write(" = ")
                                if el[0]["parent"] and val["child"]:
                                    if "str" in var_type:
                                        file.write(val["child"])
                                    elif "None" in var_type and (val["child"] is None or val["child"] == "None"):
                                        file.write(val["child"])
                                    else:
                                        if ";" in val["child"]:
                                            file.write(str(val["child"].split(";")))
                                        else:
                                            file.write(str([val["child"]]))
                                else:
                                    file.write("None")
                                file.write("\n")
                        # List handling
                        else:
                            if flag:
                                var_type = flags_dict[step][param]["Type"]
                            else:
                                var_type = params_dict[current_flag][current_flag_value][param]["Type"]
                            if all(v is None for v in el):
                                file.write("None")
                            elif "str" in var_type:
                                for ind,sel in enumerate(el):
                                    file.write(str(sel))
                                    if ind<len(el)-1:
                                        file.write(";")
                            elif "tuple" in var_type:
                                if el[0] is None or el[0] == "-inf" or el[0] == -np.inf:
                                    print(el[0])
                                    el[0] = '-2e308'
                                elif el[1] is None or el[1] == "inf" or el[1] == np.inf:
                                    el[1] = '2e308'
                                file.write(str(tuple(el)).replace("'",""))
                            else:
                                file.write(str(el))
                            file.write("\n")
    return

def write_output(workflow_data, workdir, flags_dict, params_dict, rerun=False, rerun_start=None, rerun_counter=0, docker_setup=False):
    """
    Writes all data to disk.
    
    Parameters
    ----------
    workflow_data : list
        Workflow setup.
    workdir : str
        Working directory.
    flags_dict : dict
        Dictionary for flag variables.
    params_dict : dict
        Dictionary for parameters.
    rerun : bool
        If True, the configuration is written for a rerun. The default is False
    rerun_start : str
        Step of prior run from which rerun should be started.
    rerun_counter : int
        Last counter of the prior run. The default is 0.
    docker_setup : bool
        If True, pathes are adjusted to match NextFlow setup. The default is False.
        
    Returns
    -------
    None
    
    """
    steps, precursor_metadata, language, methods, step_collection = get_steps(workflow_data, params_dict, rerun_counter)
    step_names = write_config(workdir, steps, language, rerun, rerun_start)
    write_evaluation_config(workdir, step_names, methods, rerun=rerun, rerun_start=rerun_start, docker_setup=docker_setup)
    write_parameter_file(workdir, workflow_data, step_names, step_collection, methods, precursor_metadata, flags_dict, params_dict)
    return

