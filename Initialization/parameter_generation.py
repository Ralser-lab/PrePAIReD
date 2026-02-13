#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 16:15:58 2025

@author: oliverlemke
"""

import argparse
import os
import re
import sys
import pickle
from tkinter.filedialog import askopenfilename, askdirectory

# For internal use or updates in parameter file change to False
load_dictionaries = True ####
write_to_disc = False

if not load_dictionaries:
    import pandas as pd
    import numpy as np
else:
    write_to_disc = False

#%%

# Parser to access from terminal
def init_module():
    parser = argparse.ArgumentParser(description="Step- and Method-dependent parameter setup")
    parser.add_argument("-f","--input_order", help="File containing order of steps seperated by (, or ; or space or tab or newline)",required=True)
    parser.add_argument("-i","--dict_path", help="Path of the parameter dictionaries",required=False,default=os.path.split(os.path.abspath(__file__))[0])
    parser.add_argument("-w","--workdir", help="Path to working directory where results are saved",required=False, default=".")
    parser.add_argument("-p","--output_file", help="Name of the Output parameter file",required=False, default="Parameters.txt")
    parser.add_argument("-c","--counter", help="Counter of the first step",required=False, default="1")
    parser.add_argument("-t","--parameter_overview_table", help="Source containing all parameter information",required=False, default="Parameter_Overview.tsv")
    args = parser.parse_args()
    
    # Get input
    try:
        with open(args.input_order,"r") as file:
            setup_str = file.read()
    except:
        raise ValueError("###ERROR###: Input file not found")
        sys.exit(1)
    # Set output
    outpath = args.workdir
        
    if not os.path.isdir(outpath):
        raise ValueError("###ERROR###: Output path not existing")  
        sys.exit(1)
    
    outfile = os.path.join(outpath,args.output_file)
    dict_path = args.dict_path
    config_file_path = os.path.join(outpath,"Config.txt")
    config_evaluation_file_path = os.path.join(outpath,"Config_evaluation.txt")
    parameter_overview_table_file_path = args.parameter_overview_table
    
    try:
        counter = int(args.counter)
    except:
        raise TypeError("###ERROR###: Counter need to be an integer")
    
    setup_order = get_setup_order(setup_str, counter)
    
    return setup_order, outfile, config_file_path, config_evaluation_file_path, dict_path, parameter_overview_table_file_path, outpath, counter

def get_setup_order(setup_str, counter):
    # Account for different seperators
    setup_str = setup_str.replace(" ",",").replace("\t",",").replace(";",",").replace("\n",",").replace("/",",").replace("\\",",")
    while ",," in setup_str:
        setup_str = setup_str.replace(",,",",")
    setup_str = setup_str.replace("NORMALIZATION","NORMALISATION").replace("SUMMARIZATION","SUMMARISATION").replace("FILTERINGLONG","FILTERING_LONG").replace("FILTERINGWIDE","FILTERING_WIDE")
    # Get step order
    setup_order = setup_str.split(",")
    # Remove empty lines
    setup_order = ["{:02.0f}_".format(int(counter)+ind)+item for ind,item in enumerate(setup_order) if item]
    # Add general parameters -> Will be added to .json directly later
    setup_order.insert(0,"00_GENERAL")
    return setup_order

#%%

def load_data(parameter_overview_table_file_path):
    # Read in parameter overview
    data = pd.read_csv(parameter_overview_table_file_path,sep="\t")
    # Get flag parameters
    data_flag = data[data["IsFlag"]==1]
    # Get non-flag parameters
    data_parameters = data[data["IsFlag"]==0]
    return data, data_flag, data_parameters

def get_dictionary_flags(data_flag):
    dict_flags = {step:{
        data_flag.loc[index,"Variable"]:{"Options":[el if re.search(";",str(item)) else item
                                                    for item in data_flag[(data_flag["Variable"]==data_flag.loc[index,"Variable"]) & (data_flag.loc[index,"Step"]==step)]["AllowedValues"].values 
                                                        for el in str(item).replace(" ","").split(";")],
                                         "Default":data_flag.loc[index,"PrimaryDefault"],
                                         "Comment":data_flag.loc[index,"Comment in parameters.txt"],
                                         "Type":data_flag.loc[index,"Type"],
                                         "Check Operation":data_flag.loc[index,"Operations"]
                                         } if not data_flag.loc[index,"PrimaryDefault"] is np.nan else {
                                            "Options":[el if re.search(";",str(item)) else item
                                                       for item in data_flag[(data_flag["Variable"]==data_flag.loc[index,"Variable"]) & (data_flag.loc[index,"Step"]==step)]["AllowedValues"].values 
                                                           for el in str(item).replace(" ","").split(";")],
                                            "Default":None,
                                            "Comment":data_flag.loc[index,"Comment in parameters.txt"],
                                            "Type":data_flag.loc[index,"Type"],
                                            "Check Operation":data_flag.loc[index,"Operations"]
                                            }  
                                         
        for index in data_flag[(data_flag["Step"]==step)].index} 
        for step in data_flag["Step"]}
    
    # Setup imputation setup to separate one vs. type-specific method
    dict_flags["IMPUTATION_setup"] = {"simple_imputation":{"Options":["1","0"],
                                                           "Default":"1",
                                                           "Comment":"use one (1) or missingness type specific method (0)"},
                                      "imputation_method":{"Options":dict_flags["IMPUTATION"]["imputation_method"]["Options"],
                                                           "Default":dict_flags["IMPUTATION"]["imputation_method"]["Default"],
                                                           "Comment": "simple_imputation == 0 only implemented for "+str(dict_flags["IMPUTATION"]["imputation_method_MAR"]["Options"])}
                                      }
    return dict_flags

def get_dictionary_parameters(data_parameters, data_flag):
    dict_parameters = {flag:{flag_value:{parameter:{"Options":[el if re.search(";",str(item)) else item
                                                for item in data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)]["AllowedValues"].values 
                                                    for el in str(item).replace(" ","").split(";")],
                                                    "Default":data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"PrimaryDefault"], ########
                                                    "Comment":data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"Comment in parameters.txt"],
                                                    "Type":data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"Type"],
                                                    "Check Operation":data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"Operations"]
                                                    } if not data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"PrimaryDefault"] is np.nan else {
                                                        "Options":[el if re.search(";",str(item)) else item
                                                    for item in data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)]["AllowedValues"].values 
                                                        for el in str(item).replace(" ","").split(";")],
                                                        "Default":None, ########
                                                        "Comment":data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"Comment in parameters.txt"],
                                                        "Type":data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"Type"],
                                                        "Check Operation":data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"Operations"]
                                                        }
                                for parameter in data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)]["Variable"].values}
                            for flag_value in set(data_parameters[data_parameters["Flag_variable"]==flag]["Flag_value"].values)}
                       for flag in data_flag["Variable"].values}
    return dict_parameters

def get_language(dict_parameters, data_flag):
    for flag in dict_parameters:
        for flag_value in dict_parameters[flag]:
            try:
                dict_parameters[flag][flag_value].update({"Language":data_flag.loc[data_flag[(data_flag["Variable"]==flag)&(data_flag["AllowedValues"]==flag_value)].index[0],"Language"].strip()})
            except:
                pass
    # Check if entries are missing
    for flag in dict_parameters:
        values = [value for value in data_flag[data_flag["Variable"]==flag]["AllowedValues"].values if ";" not in value]
        try:
            for flag_value in values:
                if not flag_value in dict_parameters[flag]:
                    dict_parameters[flag][flag_value]={"Language":data_flag.loc[data_flag[(data_flag["Variable"]==flag)&(data_flag["AllowedValues"]==flag_value)].index[0],"Language"].strip()}
        except:
            print(flag_value)
    
    
    return dict_parameters
    
def expand_data_parameters(data_parameters):
    # Expand data_parameters to break down ";" in Flag_value/variable
    indices = data_parameters.index
    for index in indices:
        if re.search(";",str(data_parameters.loc[index,"Flag_value"])):
            values = data_parameters.loc[index,"Flag_value"].replace(" ","").split(";")
            for value in values:
                new_line = pd.DataFrame(data_parameters.loc[index].copy()).T
                new_line.index = ["tmp"]
                new_line["Flag_value"] = value
                data_parameters = pd.concat([data_parameters,new_line],axis=0)
    data_parameters.reset_index(drop=True,inplace=True)
    
    indices = data_parameters.index
    for index in indices:
        if re.search(";",str(data_parameters.loc[index,"Flag_variable"])):
            values = data_parameters.loc[index,"Flag_variable"].replace(" ","").split(";")
            for value in values:
                new_line = pd.DataFrame(data_parameters.loc[index].copy()).T
                new_line.index = ["tmp"]
                new_line["Flag_variable"] = value
                data_parameters = pd.concat([data_parameters,new_line],axis=0)
    data_parameters.reset_index(drop=True,inplace=True)

    # Drop expanded rows
    data_parameters.drop(index=[index for index in data_parameters.index 
                                  if re.search(";",str(data_parameters.loc[index,"Flag_value"]))
                                                       or re.search(";",str(data_parameters.loc[index,"Flag_variable"]))],inplace=True)
    data_parameters.reset_index(drop=True,inplace=True)
    return data_parameters
    
def expand_data_flags(data_flag):
    indices = data_flag.index
    for index in indices:
        if re.search(";",str(data_flag.loc[index,"AllowedValues"])):
            values = data_flag.loc[index,"AllowedValues"].replace(" ","").split(";")
            for value in values:
                new_line = pd.DataFrame(data_flag.loc[index].copy()).T
                new_line.index = ["tmp"]
                new_line["AllowedValues"] = value
                data_flag = pd.concat([data_flag,new_line],axis=0)
    data_flag.reset_index(drop=True,inplace=True)
    return data_flag

def generate_dictionaries(parameter_overview_table_file_path):
    # Read in parameter overview
    data, data_flag, data_parameters = load_data(parameter_overview_table_file_path)

    # Get dictionary flags (step:flag:flag_values->options,default,comment)
    dict_flags =  get_dictionary_flags(data_flag)

    # Expand data_parameters to break down ";" in Flag_value/variable
    data_parameters = expand_data_parameters(data_parameters)

    # Get dictionary parameters (flag:flag_value:parameter->options,default,comment)
    dict_parameters = get_dictionary_parameters(data_parameters, data_flag)
    
    # Expand data_flags with all allowed flag values (needed for language)
    data_flag = expand_data_flags(data_flag)

    # Get programming language property
    dict_parameters = get_language(dict_parameters, data_flag)
    
    return dict_flags, dict_parameters

#%%
def check_steps(setup_order, dict_flags):
    # Check if all setup steps are known
    check = [el[3:] for el in setup_order if el[3:] not in dict_flags]
    if check:
        raise ValueError("###ERROR###: Preprocessing step "+str(check)+" not found.")
        sys.exit(1)
        
def check_language(language):
    for item in language:
        if len(set(item))>1:
            raise ValueError("###Internal Error###: Please report and setup "+config_file_path+" manually according to documentation.")
            sys.exit(1)
            
#%%
        
def write_flag_parameter_file(setup_order, outfile, dict_flags):
    # Save to file with Imputation replaced by Imputation_setup
    with open(outfile,"w") as file:
        for setup_c in setup_order:
            setup = setup_c[3:]
            if setup == "IMPUTATION":
                setup = "IMPUTATION_setup"
                file.write("%%%%%%%%%%"+re.search("[0-9]+",setup_c).group(0)+"_"+setup+"%%%%%%%%%%\n\n")
            else:
                file.write("%%%%%%%%%%"+setup_c+"%%%%%%%%%%\n\n")
            for flag in dict_flags[setup]:
                if type(dict_flags[setup][flag]["Comment"]) is str:
                    file.write("# "+dict_flags[setup][flag]["Comment"]+"\n")
                #file.write("## Options: "+str(dict_flags[setup][flag]["Options"]).replace("\'","")+"\n")          
                #if dict_flags[setup][flag]["Options"] == [np.nan]:
                #    pass
                #elif
                if dict_flags[setup][flag]["Options"][0].startswith("(") or dict_flags[setup][flag]["Options"][0].startswith("["):
                    file.write("## Options: range"+str(dict_flags[setup][flag]["Options"][0]).replace("\'","")+"\n") 
                elif len(dict_flags[setup][flag]["Options"]) == 1:
                    file.write("## Options: "+str(dict_flags[setup][flag]["Options"][0]).replace("\'","")+"\n") 
                else:
                    file.write("## Options: "+str(dict_flags[setup][flag]["Options"]).replace("\'","")+"\n") 
                file.write(flag+" = "+str(dict_flags[setup][flag]["Default"])+"\n\n")
            file.write("\n")
            
#%%

def get_blocks(outfile):
    with open(outfile,"r") as file:
        blocks = []
        block = []
        for ind,line in enumerate(file.read().splitlines()):
            if ind == 0:
                block.append(line)
            elif line.startswith("%%%%"):
                blocks.append(block)
                block = [line]
            else:
                block.append(line)
        blocks.append(block)
    return blocks

def write_parameter_file_imputation(outfile, dict_flags):
    # Read in and extract blocks
    blocks = get_blocks(outfile)
    
    # Replace the IMPUTATION_setup block with IMPUTATION Block and flags (one vs. type-specific method)
    for ind,block in enumerate(blocks):
        if "IMPUTATION_setup%%%%%%%%%%" in block[0]:
            block_new = ["%%%%%%%%%%"+re.search("[0-9]+",block[0]).group(0)+"_IMPUTATION%%%%%%%%%%"]
            block_new.append("")
            for line in block:
                # check imputation type
                if line.startswith("simple_imputation"):
                    if line.split(" = ")[1].strip()=="1":
                        simple = True
                    else:
                        simple = False
                    block_new.append("## FLAG DO NOT TOUCH ")
                    block_new.append(line)
            # extract imputation method and check for existence
            for line in block:
                if line.startswith("imputation_method"):
                    imputation_method = line.split(" = ")[1].strip()
                    if not simple:
                        if imputation_method not in dict_flags["IMPUTATION"]["imputation_method_MAR"]["Options"]:
                            raise ValueError("###ERROR###: Method not implemented for complex setup")
                            sys.exit(1) 
            setup = "IMPUTATION"
            block_new.append("## FLAG DO NOT TOUCH ")
            block_new.append("imputation_method = "+imputation_method)
            block_new.append("")
            for flag in dict_flags[setup]:
                if simple:
                    if imputation_method in dict_flags["IMPUTATION"]["imputation_method"]["Options"]:
                        # Distinguish between R and Python: Python-specific flags
                        if imputation_method not in ["missforest","bpca_R","ppca","ppm","MNAR"]:
                            if flag not in ["imputation_method_MAR","imputation_method_MNAR","imputation_method_ProtAbsence","imputation_method_ProtAbsenceMAR","imputation_method_ProtAbsenceMNAR","simple_imputation","imputation_method"]:
                                if type(dict_flags[setup][flag]["Comment"]) is str: 
                                    block_new.append("# "+str(dict_flags[setup][flag]["Comment"]))
                                block_new.append("## Options: "+str(dict_flags[setup][flag]["Options"]).replace("\'",""))
                                block_new.append(flag+" = "+str(dict_flags[setup][flag]["Default"])+"\n")
                        # R
                        elif flag in ["imputation_method"]:
                            block_new.append("## Options: "+str(dict_flags["IMPUTATION"]["imputation_method"]["Options"]).replace("\'",""))
                            block_new.append("imputation_method = "+imputation_method+"\n")
                    else:
                        raise ValueError("###ERROR###: Method for IMPUTATION not found; Allowed values in "+str(dict_flags["IMPUTATION"]["imputation_method"]["Options"]+dict_flags["IMPUTATION"]["imputation_method_MAR"]["Options"]).replace("\'",""))
                        sys.exit(1)             
                else:
                    # Setup type specific flags (Python)
                    if flag not in ["imputation_method","simple_imputation"]:
                        if type(dict_flags[setup][flag]["Comment"]) is str: 
                            block_new.append("# "+str(dict_flags[setup][flag]["Comment"]))    
                        if dict_flags[setup][flag]["Options"][0].startswith("(") or dict_flags[setup][flag]["Options"][0].startswith("["):
                            block_new.append("## Options: range"+str(dict_flags[setup][flag]["Options"][0]).replace("\'","")) 
                        elif len(dict_flags[setup][flag]["Options"]) == 1:
                            block_new.append("## Options: "+str(dict_flags[setup][flag]["Options"][0]).replace("\'","")) 
                        else:
                            block_new.append("## Options: "+str(dict_flags[setup][flag]["Options"]).replace("\'","")) 
                        block_new.append(flag+" = "+str(dict_flags[setup][flag]["Default"])+"\n") 
            block_new.append("")
            blocks[ind] = block_new
    # Save to file
    write_parameter_file(outfile,blocks)

def write_parameter_file(outfile,blocks):
    with open(outfile,"w") as file:
        for ind,block in enumerate(blocks): 
            for ind_line,line in enumerate(block):
                # Avoid empty lines at the end of the parameter file
                if ind == len(blocks)-1 and ind_line  == len(block)-1:
                    file.write(line)
                else:
                    file.write(line+"\n")

#%%                    
def write_final_parameter_file(outfile,dict_flags,dict_parameters):
    # Read in and extract blocks
    blocks = get_blocks(outfile) 
    
    # Save language per module
    language = []
    # Add parameters for method to the blocks
    for ind,block in enumerate(blocks): 
        # Reset block
        block_new = []
        # Avoid duplicated parameters in parameter file
        parameters_included = []
        language_coll = []
        for ind_line,line in enumerate(block):
            # Extract step
            if ind_line == 0:
                setup = re.match("%%%%%%%%%%(.+?)%%%%%%%%%%",block[0]).group(1)[3:]   
            # Add comment so that flags are not touched anymore
            if line.startswith("## Options"):
                block_new.append("## FLAG DO NOT TOUCH ")
            else:
                block_new.append(line)
            # get flags and flag_values
            if line and not line.startswith("%%%") and not line.startswith("#"):
                flag, v_ = line.split(" = ")
                # ensure some flags can have several flag_values
                if line.startswith("filtering_wide_methods") or line.startswith("imputation_sample_hierarchies"):
                    flag_values = v_.split(";")
                else:
                    flag_values = [v_]
                
                # Get parameters for flag:flag_value combinations
                for flag_value in flag_values:
                    # Check whether flag value exists
                    if setup != "GENERAL" and not flag_value in dict_flags[setup][flag]["Options"] and not any([re.search("column",el) for el in dict_flags[setup][flag]["Options"]]):
                        raise ValueError("###ERROR###: "+str(flag_value) +" for "+str(flag)+" not found.")
                        sys.exit(1)
                    
                    # Avoid errors due to missing dict entries for standalone flags
                    if dict_parameters[flag]:
                        if flag_value in dict_parameters[flag]:
                            for parameter,param_item in dict_parameters[flag][flag_value].items():
                                # Skip language imformation for parameter file (needed later)
                                if param_item not in ["R","Python"]:
                                    # Check for redundancies in module e.g. qvalue_treshold
                                    if parameter not in parameters_included:
                                        if type(param_item["Comment"]) is str: 
                                            block_new.append("# "+str(param_item["Comment"]))
                                        #if param_item["Options"] == [np.nan]:
                                        #    pass
                                        #elif 
                                        if param_item["Options"][0].startswith("(") or param_item["Options"][0].startswith("["):
                                            block_new.append("## Options: range"+str(param_item["Options"][0]).replace("\'","")) 
                                        elif len(param_item["Options"]) == 1:
                                            block_new.append("## Options: "+str(param_item["Options"][0]).replace("\'","")) 
                                        else:
                                            block_new.append("## Options: "+str(param_item["Options"]).replace("\'","")) 
                                        block_new.append(str(parameter)+" = "+str(param_item["Default"]))
                                        parameters_included.append(parameter)
                                else:
                                    language_coll.append(param_item)
        # Add block
        blocks[ind] = block_new
        # Add language collection
        language.append(language_coll)

    # Save to file
    write_parameter_file(outfile,blocks)   
    
    return language

def write_config(config_file_path, setup_order, language):
    with open(config_file_path,"w") as file:
        file.write("Run = new\n")
        file.write("order = ")
        c = 0
        for step, language_step in zip(setup_order,language):
            if step != "00_GENERAL":
                if "IMPUTATION" in step:
                    if language_step[1] == "Python":
                        file.write(step+"_py")
                    elif language_step[1] == "R":
                        file.write(step+"_R")
                    else:
                        raise ValueError("###Internal Error###: Please report and fix Config.txt manually according to documentation")
                        sys.exit(1)
                else:
                    if language_step[0] == "Python":
                        file.write(step+"_py")
                    elif language_step[0] == "R":
                        file.write(step+"_R")
                    else:
                        raise ValueError("###Internal Error###: Please report and fix Config.txt manually according to documentation")
                        sys.exit(1)
            c+=1
            if c<len(setup_order) and step != "00_GENERAL":
                file.write(", ")

def write_evaluation_config(config_evaluation_file_path, setup_order,blocks_flags, outpath, counter):
    blocks_config_eval = []
    for ind,block in enumerate(blocks_flags):
        counter_ = int(re.search("\%\%([0-9]+?)\_",block[0]).group(1))
        step = re.search("\%\%[0-9][0-9]\_(.+?)\%",block[0]).group(1)
        block_new = [block[0]]
        block_new.append("%-%-%-%")
        if counter_ == 0:
            for line in block:
                if line and not line.startswith("#") and not line.startswith("%"):
                    block_new.append(line)
            block_new.append("param_fp = "+os.path.join(outpath,"Parameter_Files","00_GENERAL_Parameters.json"))
            blocks_config_eval.append(block_new)   
        else:
            #block_new = [block[0]]
            ##### Check if layout works
            if ind > 1:
                block_new.append("data_origin = "+os.path.join(outpath,setup_order[ind-1],"data_"+"{:02d}".format(ind-1+counter-1)+".tsv"))
                #block_new.append("mask_description_origin = "+os.path.join(outpath,setup_order[ind-1],"data_"+"{:02d}".format(ind-1)+"_mask_desription.tsv"))
                #block_new.append("data_processed = "+os.path.join(outpath,setup_order[ind],"data_"+"{:02d}".format(ind)+".tsv"))
                #block_new.append("mask_processed = "+os.path.join(outpath,setup_order[ind],"data_"+"{:02d}".format(ind)+"_mask.tsv"))
                #block_new.append("mask_description_processed = "+os.path.join(outpath,setup_order[ind],"data_"+"{:02d}".format(ind)+"_mask_desription.tsv"))
            else:
                block_new.append("data_origin = None")
            block_new.append("data_processed = "+os.path.join(outpath,setup_order[ind],"data_"+"{:02d}".format(ind+counter-1)+".tsv"))
                #block_new.append("mask_description_origin = None")
                #block_new.append("mask_processed = "+os.path.join(outpath,setup_order[ind],"data_"+"{:02d}".format(ind)+"_mask.tsv"))
                #block_new.append("mask_description_processed = "+os.path.join(outpath,setup_order[ind],"data_"+"{:02d}".format(ind)+"_mask_desription.tsv"))
            block_new.append(("param_fp = "+os.path.join(outpath,"Parameter_Files",setup_order[ind]+"_Parameters.json")))
            block_new.append("%-%-%-%")
            #
            if step in ["FILTERING_LONG","FILTERING_WIDE","NORMALISATION","IMPUTATION","CORRECTION","SUMMARISATION"]:
                block_new.append("plot_int_box = 1")
                ### Use F notation for figure legend
                block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F01.png"))
                if step == "FILTERING_LONG":
                    block_new.append("size = (5,10)")
                else:
                    block_new.append("size = (10,10)")
                ### Coloring column needed
                #if step == "CORRECTION":
                #    block_new.append("use_metadata_coloring = 1")
                #    block_new.append("use_parameter_file = 1")
            else:
                block_new.append("plot_int_box = 0")
                block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F01.png"))
            block_new.append("%-%-%-%")
            #  
            if step in ["FILTERING_LONG","FILTERING_WIDE"]:
                block_new.append("plot_filter_mask = 1")
            else:
                block_new.append("plot_filter_mask = 0")
            block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F02.png"))
            block_new.append("%-%-%-%")
            #
            if step in ["NORMALISATION","IMPUTATION","CORRECTION","SUMMARISATION"]:
                block_new.append("plot_cv = 1")
            else:
                block_new.append("plot_cv = 0")
            block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F03.png"))
            block_new.append("%-%-%-%")
            ### column_color from parameter_file
            ### for correction decide on property if it is continuous or not -> continuous = None
            ### for summarization -> sample_type
            if step in ["NORMALISATION","IMPUTATION","CORRECTION","SUMMARISATION"]:
                block_new.append("plot_pca = 1")
                block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F04.png"))
                if step == "NORMALISATION":
                    block_new.append("continuous = True")
                    block_new.append("coloring_property = Missingness")
                    block_new.append("use_metadata_coloring = False")
                    #block_new.append("use_parameter_file = 0")
                elif step == "IMPUTATION":
                    block_new.append("continuous = True")
                    block_new.append("coloring_property = Mean_Intensities")
                    block_new.append("use_metadata_coloring = False")
                    #block_new.append("use_parameter_file = 0")
                elif step == "CORRECTION":
                    block_new.append("continuous = None")
                    block_new.append("coloring_property = None")
                    block_new.append("use_metadata_coloring = True")
                    #block_new.append("use_parameter_file = 1")
                elif step == "SUMMARISATION":
                    block_new.append("continuous = False")
                    block_new.append("coloring_property = None")
                    block_new.append("use_metadata_coloring = True")
                    #block_new.append("use_parameter_file = 1")  
            else:
                block_new.append("plot_pca = 0")
                block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F04.png"))
                block_new.append("continuous = True")
                block_new.append("coloring_property = Missingness")
                block_new.append("use_metadata_coloring = False")
            block_new.append("%-%-%-%")
            #    
            if step in ["NORMALISATION"]:
                block_new.append("plot_distribution = 1")
            else:
                block_new.append("plot_distribution = 0")
            block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F05.png"))
            block_new.append("%-%-%-%")
            ### Confounder (first entry), for batchwise color according to column and remove linreg (or adjust linreg)
            if step in ["NORMALISATION","CORRECTION"]:
                if step == "NORMALISATION":
                    block_new.append("plot_scatter_series = 1")
                    #block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F06.png"))
                    #block_new.append("plot_by = Correlation")
                    #block_new.append("use_parameter_file = 0")
                else:
                    c = False
                    for line in block:
                        if "correct_method = adjust_confounders" in line or "correct_method = adjust_confounders_batchwise" in line:
                            block_new.append("plot_scatter_series = 1")
                            #block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F06.png"))
                            #block_new.append("plot_by = Correlation")
                            #block_new.append("use_parameter_file = 1")
                            c = True
                    if not c:
                        block_new.append("plot_scatter_series = 0")
            else:
                block_new.append("plot_scatter_series = 0")     
            block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F06.png"))
            block_new.append("plot_by = Correlation")
            block_new.append("%-%-%-%")
            #### Method 
            #### Cutoff (intensity_dependent_missingness_estimation) 
            if step in ["IMPUTATION"]: 
                c = False
                for line in block:
                    if line == "eval_missingness_by_detection_threshold = 1":
                        block_new.append("plot_decision_boundary = 1")
                        #block_new.append("use_parameter_file = 1")
                        c = True
                if not c:
                    block_new.append("plot_decision_boundary = 0")
            else:
                block_new.append("plot_decision_boundary = 0")
            block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F07.png"))
            block_new.append("%-%-%-%")
            #### Correlation type
            #### Confounder (first entry)
            if step in ["CORRECTION"]:
                c = False
                for line in block:
                    if "correct_method = adjust_confounders" in line or "correct_method = adjust_confounders_batchwise" in line:
                        block_new.append("plot_correlation = 1")
                        #block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F08.png"))
                        #block_new.append("use_parameter_file = 1")
                        c = True
                if not c:
                    block_new.append("plot_correlation = 0")
            else:
                block_new.append("plot_correlation = 0")
            block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F08.png"))
            block_new.append("%-%-%-%")
            #
            if step in ["NORMALISATION"]: 
                c = False
                for line in block:
                    if "normalization_method = median_scaling" in line:
                        block_new.append("plot_distribution_hist = 1")
                        block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F09.png"))
                        block_new.append("normalization_factors = "+os.path.join(outpath,setup_order[ind],"Normalization_factors.tsv"))
                        c = True
                if not c:
                    block_new.append("plot_distribution_hist = 0")
                    block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F09.png"))
            else:
                block_new.append("plot_distribution_hist = 0")
                block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","report_folder","F09.png"))
            block_new.append("%-%-%-%")
            #
            if step in ["FILTERING_LONG","FILTERING_WIDE"]: 
                block_new.append("plot_abundance_vs_completeness = 1")
                block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F10.png"))
                if step == "FILTERING_LONG":
                    block_new.append("gridsize = 50")
                    block_new.append("n_bins = 50")
                elif step == "FILTERING_WIDE":
                    block_new.append("gridsize = 25")
                    block_new.append("n_bins = 25")
            else:
                block_new.append("plot_abundance_vs_completeness = 0") 
                block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F10.png"))
            block_new.append("%-%-%-%")
            #
            if step in ["FILTERING_LONG","FILTERING_WIDE"]: 
                block_new.append("plot_overall_completeness = 1")
                #block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F11.png"))
            else:
                block_new.append("plot_overall_completeness = 0") 
            block_new.append("outfile = "+ os.path.join(outpath,setup_order[ind],"report_folder","F11.png"))
            #"plot_imputation_vs_covariate = 0" ### Skip for now as this needs parameters from multiple steps
            #block_new.append("%-%-%-%")
            blocks_config_eval.append(block_new)   
            
    with open(config_evaluation_file_path,"w") as file:
        for block in blocks_config_eval:
            for line in block:
                file.write(line+"\n")
    return
        
#%%

if __name__ == "__main__":
    
    # Get input arguments
    if len(sys.argv) <= 1:
        
        counter = int(input("Please define the counter value of the first step (int): "))
        if counter <= 0:
            raise ValueError("Only counters bigger than 0 are allowed.")
        
        input_order = askopenfilename(title="Please select the file containing the setup order")
        with open(input_order,"r") as file:
            setup_str = file.read()
        
        outpath = askdirectory(title="Please select your working directory", initialdir=os.path.abspath(os.getcwd()))
        outfile = os.path.join(outpath,"Parameters.txt")
        config_file_path = os.path.join(outpath,"Config.txt")
        config_evaluation_file_path = os.path.join(outpath,"Config_evaluation.txt")
        
        if load_dictionaries:
            dict_path = askdirectory(title="Please select the dictionary path", initialdir=os.path.abspath(os.getcwd()))
        else:
            parameter_overview_table_file_path = askopenfilename(title="Please select the parameter overview file")
            if write_to_disc:
                dict_path = askdirectory(title="Please select the path for saving the dictionaries", initialdir=os.path.abspath(os.getcwd()))
            else:
                dict_path = None
        
        setup_order = get_setup_order(setup_str, counter)
    
    else:
        setup_order, outfile, config_file_path, config_evaluation_file_path, dict_path, parameter_overview_table_file_path, outpath, counter = init_module()
    
    # Get dictionaries
    if not load_dictionaries:
        dict_flags, dict_parameters = generate_dictionaries(parameter_overview_table_file_path)
        if write_to_disc:
            pickle.dump(dict_flags,open(os.path.join(dict_path,"dict_flags.pkl"),"wb"))
            pickle.dump(dict_parameters,open(os.path.join(dict_path,"dict_parameters.pkl"),"wb"))
    else:
        dict_flags = pickle.load(open(os.path.join(dict_path,"dict_flags.pkl"),"rb"))
        dict_parameters = pickle.load(open(os.path.join(dict_path,"dict_parameters.pkl"),"rb"))
    
    # Check for internal error
    check_steps(setup_order, dict_flags)
    
    write_flag_parameter_file(setup_order, outfile, dict_flags)
    # Raw input -> User can check flags and imputation method
    input("Please set the methods for every step. The parameters file is saved in "+outfile+".\n\tPress Enter after the methods are set.")
    
    # Save to file
    write_parameter_file_imputation(outfile, dict_flags)
    # Raw input -> User can check flags
    input("Please set the imputation method and flags if necessary. If no imputation step is used, ignore this step. The parameters file is saved in "+outfile+".\n\tPress Enter after the methods are set.")
                  
    # Save to file
    blocks_flags = get_blocks(outfile)
    language = write_final_parameter_file(outfile,dict_flags,dict_parameters)
    print("The parameters file is saved in "+outfile+".")
    
    # Check for internal error
    #check_language(language)
    
    # Create config file
    write_config(config_file_path, setup_order, language)
    print("The configuration file is saved in "+config_file_path+".")
    
    # Create evaluation config file
    write_evaluation_config(config_evaluation_file_path, setup_order, blocks_flags, outpath, counter)
    print("The configuration file for evaluation is saved in "+config_file_path+".")
