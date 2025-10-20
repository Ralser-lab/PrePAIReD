#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 16:15:58 2025

@author: oliverlemke
"""

import pandas as pd
import argparse
import os
import re
import sys

#%%
###### TODO ######
### Setup config -> DONE
##### Check of parser and config need to deal with run/rerun -> DONE
### Write function to check if all parameters are present -> Implement in parameter check
### Add parser -> DONE

### Add value if no imputation should be done like "no_imputation" -> Add to table
#%%

# Parser to access from terminal
#"""
parser = argparse.ArgumentParser(description="Step- and Method-dependent parameter setup")
parser.add_argument("-f","--input_order", help="File containing order of steps seperated by (, or ; or space or tab or newline)",required=True)
parser.add_argument("-o","--output_path", help="Path for the Output parameter file",required=False, default=".")
args = parser.parse_args()

# Get input
try:
    with open(args.input_order,"r") as file:
        setup_str = file.read()
except:
    raise ValueError("###ERROR###: Input file not found")
    sys.exit(1)
# Set output‚
outpath = args.output_path
    
if not os.path.isdir(outpath):
    raise ValueError("###ERROR###: Output path not existing")  
    sys.exit(1)

outfile = os.path.join(outpath,"Parameters_test.txt")
config_file_path = os.path.join(outpath,"Config.txt")
#"""
"""
setup_str = "FILTERING_LONG,FILTERING_WIDE,NORMALIZATION,IMPUTATION,CORRECTION,NORMALIZATION,SUMMARIZATION"
outfile = "Parameters_test.txt"
config_file_path = "Config.txt"
"""
# Account for different seperators
setup_str = setup_str.replace(" ",",").replace("\t",",").replace(";",",").replace("\n",",").replace(",,",",").replace(",,",",")
# Get step order
setup_order = setup_str.split(",")
# Remove empty lines
setup_order = [item for item in setup_order if item]

#%%
# Add general parameters -> Will be added to .json directly later
setup_order.insert(0,"GENERAL")

# Read in parameter overview
data = pd.read_csv("Preprocessing Working Tasks - Parameter Overview.tsv",sep="\t")
# Get flag parameters
data_flag = data[data["IsFlag"]==1]
# Get non-flag parameters
data_parameters = data[data["IsFlag"]==0]

# Check if all setup steps are known
check = [el for el in setup_order if el not in data_flag["Step"].values]
if check:
    raise ValueError("###ERROR###: Preprocessing step not found.")
    sys.exit(1)

#%%
# Get dictionary flags (step:flag:flag_values->options,default,comment)
dict_flags = {step:{
    data_flag.loc[index,"Variable"]:{"Options":[el if re.search(";",str(item)) else item
                                                for item in data_flag[(data_flag["Variable"]==data_flag.loc[index,"Variable"]) & (data_flag.loc[index,"Step"]==step)]["AllowedValues"].values 
                                                    for el in str(item).replace(" ","").split(";")],
                                     "Default":data_flag.loc[index,"PrimaryDefault"],
                                     "Comment":data_flag.loc[index,"Comment in parameters.txt"]
                                     }
    for index in data_flag[(data_flag["Step"]==step)].index} 
    for step in data_flag["Step"]}

# Setup imputation setup to separate one vs. type-specific method
dict_flags["IMPUTATION_setup"] = {"simple_imputation":{"Options":[1,0],
                                                            "Default":1,
                                                            "Comment":"use one (1) or missingness type specific method (0)"},
                                        "imputation_method":{"Options":dict_flags["IMPUTATION"]["imputation_method"]["Options"],
                                                             "Default":dict_flags["IMPUTATION"]["imputation_method"]["Default"],
                                                             "Comment": "imputation_setup == 0 only implemented for "+str(dict_flags["IMPUTATION"]["imputation_method_MAR"]["Options"])}
                                        }

#%%

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
#%%
# Drop expanded rows
data_parameters.drop(index=[index for index in data_parameters.index 
                              if re.search(";",str(data_parameters.loc[index,"Flag_value"]))
                                                   or re.search(";",str(data_parameters.loc[index,"Flag_variable"]))],inplace=True)
data_parameters.reset_index(drop=True,inplace=True)

#%%

# Get dictionary parameters (flag:flag_value:parameter->options,default,comment)
dict_parameters = {flag:{flag_value:{parameter:{"Options":[el if re.search(";",str(item)) else item
                                            for item in data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)]["AllowedValues"].values 
                                                for el in str(item).replace(" ","").split(";")],
                                                "Default":data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"Default"], ########
                                                "Comment":data_parameters.loc[data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)&(data_parameters["Variable"]==parameter)].index[0],"Comment in parameters.txt"]}
                            for parameter in data_parameters[(data_parameters["Flag_variable"]==flag)&(data_parameters["Flag_value"]==flag_value)]["Variable"].values}
                        for flag_value in set(data_parameters[data_parameters["Flag_variable"]==flag]["Flag_value"].values)}
                   for flag in data_flag["Variable"].values}

#%%
# Expand data_flags with all allowed flag values
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
#%%
# Get programming language property
for flag in dict_parameters:
    for flag_value in dict_parameters[flag]:
        try:
            dict_parameters[flag][flag_value].update({"Language":data_flag.loc[data_flag[(data_flag["Variable"]==flag)&(data_flag["AllowedValues"]==flag_value)].index[0],"Language"].strip()})
        except:
            pass

#%%
# Save to file with Imputation replaced by Imputation_setup
with open(outfile,"w") as file:
    for setup in setup_order:
        if setup == "IMPUTATION":
            setup = "IMPUTATION_setup"
        file.write("%%%%%%%%%%"+setup+"%%%%%%%%%%\n\n")
        for flag in dict_flags[setup]:
            if type(dict_flags[setup][flag]["Comment"]) is str:
                file.write("# "+dict_flags[setup][flag]["Comment"]+"\n")
            file.write("## Options: "+str(dict_flags[setup][flag]["Options"]).replace("\'","")+"\n")
            file.write(flag+" = "+str(dict_flags[setup][flag]["Default"])+"\n")
        file.write("\n")

#%%
# Raw input -> User can check flags and imputation method
input("Please set the methods for every step. The parameters file is saved in "+outfile+".\n\tPress Enter after the methods are set.")

#%%
# Read in and extract blocks
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

#%%

# Replace the IMPUTATION_setup block with IMPUTATION Block and flags (one vs. type-specific method)
for ind,block in enumerate(blocks):
    if block[0]=="%%%%%%%%%%IMPUTATION_setup%%%%%%%%%%":
        block_new = ["%%%%%%%%%%IMPUTATION%%%%%%%%%%"]
        block_new.append("")
        for line in block:
            # check imputation type
            if line.startswith("simple_imputation"):
                if line.split(" = ")[1]=="1":
                    simple = True
                else:
                    simple = False
                block_new.append(line)
        # extract imputation method and check for existence
        for line in block:
            if line.startswith("imputation_method"):
                imputation_method = line.split(" = ")[1]
                if not simple:
                    if imputation_method not in dict_flags["IMPUTATION"]["imputation_method_MAR"]["Options"]:
                        raise ValueError("###ERROR###: Method not implemented for complex setup")
                        sys.exit(1) 
        #if imputation_method in dict_flags["IMPUTATION"]["imputation_method"]["Options"]:
        #    if type(dict_flags[setup][flag]["Comment"]) is str: 
        #        block_new.append("# "+str(dict_flags["IMPUTATION"]["imputation_method"]["Comment"]))
        #    block_new.append("## Options: "+str(dict_flags["IMPUTATION"]["imputation_method"]["Options"]).replace("\'",""))      
        #    block_new.append("imputation_method = "+imputation_method)
        #elif imputation_method in dict_flags["IMPUTATION"]["imputation_method_MAR"]["Options"]:
        # replace Imputation_setup with Imputation
        setup = "IMPUTATION"
        block_new.append("### IMPORTANT NOTE: Imputation methods set to default. Please adjust if necessary.")
        for flag in dict_flags[setup]:
            if simple:
                if imputation_method in dict_flags["IMPUTATION"]["imputation_method"]["Options"]:
                    # Distinguish between R and Python: Python-specific flags
                    if imputation_method not in ["missforest","bpca_R","ppca","ppm","MNAR"]:
                        if flag not in ["imputation_method_MAR","imputation_method_MNAR","imputation_method_ProtAbsence","imputation_method_ProtAbsenceMAR","imputation_method_ProtAbsenceMNAR","simple_imputation"]:
                            #if flag == "imputation_method":
                            #    block_new.append("## Options: "+str(dict_flags["IMPUTATION"]["imputation_method_MAR"]["Options"]+dict_flags["IMPUTATION"]["imputation_method"]["Options"]).replace("\'",""))
                            #    block_new.append("imputation_method = "+imputation_method)
                            #else:
                            if type(dict_flags[setup][flag]["Comment"]) is str: 
                                block_new.append("# "+str(dict_flags[setup][flag]["Comment"]))
                            block_new.append("## Options: "+str(dict_flags[setup][flag]["Options"]).replace("\'",""))
                            block_new.append(flag+" = "+str(dict_flags[setup][flag]["Default"]))
                    # R
                    elif flag in ["imputation_method"]:
                        block_new.append("## Options: "+str(dict_flags["IMPUTATION"]["imputation_method"]["Options"]).replace("\'",""))
                        block_new.append("imputation_method = "+imputation_method)
                else:
                    raise ValueError("###ERROR###: Method for IMPUTATION not found; Allowed values in "+str(dict_flags["IMPUTATION"]["imputation_method"]["Options"]+dict_flags["IMPUTATION"]["imputation_method_MAR"]["Options"]).replace("\'",""))
                    sys.exit(1)             
            else:
                # Setup type specific flags (Python)
                if flag not in ["imputation_method"]:
                    #if flag in ["imputation_method_MAR","imputation_method_MNAR","imputation_method_ProtAbsence","imputation_method_ProtAbsenceMAR","imputation_method_ProtAbsenceMNAR","simple_imputation"]:
                    #    block_new.append("## Options: "+str(dict_flags[setup][flag]["Options"]).replace("\'",""))
                    #    block_new.append(flag+" = "+str(dict_flags[setup][flag]["Default"])) #### maybe change later for different default methods ##### set ProtAbsence to None or use PrimaryDefault here
                    #else:
                    if type(dict_flags[setup][flag]["Comment"]) is str: 
                        block_new.append("# "+str(dict_flags[setup][flag]["Comment"]))    
                    block_new.append("## Options: "+str(dict_flags[setup][flag]["Options"]).replace("\'","")) 
                    block_new.append(flag+" = "+str(dict_flags[setup][flag]["Default"])) 
        block_new.append("")
        blocks[ind] = block_new
      
#%%
# Save to file
with open(outfile,"w") as file:
    for ind,block in enumerate(blocks): 
        for ind_line,line in enumerate(block):
            # Avoid empty lines at the end of the parameter file
            if ind == len(blocks)-1 and ind_line  == len(block)-1:
                file.write(line)
            else:
                file.write(line+"\n")

#%%
# Raw input -> User can check flags
input("Please set the imputation method and flags if necessary. The parameters file is saved in "+outfile+".\n\tPress Enter after the methods are set.")
                    
#%%
# Read in and extract blocks
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
    
#%%
# Save language per module
language = []
# Save steps and programming language
steps = []
# Add parameters for method to the blocks
for ind,block in enumerate(blocks): 
    # Reset block
    block_new = []
    # Avoid dupolicated parameters in parameter file
    parameters_included = []
    language_coll = []
    for ind_line,line in enumerate(block):
        # Extract step
        if ind_line == 0:
            setup = re.match("%%%%%%%%%%(.+?)%%%%%%%%%%",block[0]).group(1)   
            steps.append(setup)
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
                                    block_new.append("## Options: "+str(param_item["Options"]).replace("\'","")) 
                                    block_new.append(str(parameter)+" = "+str(param_item["Default"]))
                                    parameters_included.append(parameter)
                            else:
                                language_coll.append(param_item)
    # Add block
    blocks[ind] = block_new
    # Add language collection
    language.append(language_coll)

#%%
# Save to file
with open(outfile,"w") as file:
    for ind,block in enumerate(blocks): 
        for ind_line,line in enumerate(block):
            # Avoid empty lines at the end of the parameter file
            if ind == len(blocks)-1 and ind_line  == len(block)-1:
                file.write(line)
            else:
                file.write(line+"\n")
print("The parameters file is saved in "+outfile+".")

#%%
# Check for internal error
for item in language:
    if len(set(item))>1:
        raise ValueError("###Internal Error###: Please report and setup Config.txt manually according to documentation.")
        sys.exit(1)
#%%
# Create config file
with open(config_file_path,"w") as file:
    file.write("Run = new\n")
    file.write("order = ")
    c = 0
    for step, language_step in zip(steps,language):
        if step != "GENERAL":
            if language_step[0] == "Python":
                file.write(step+"_py")
            elif language_step[0] == "R":
                file.write(step+"_R")
            else:
                raise ValueError("###Internal Error###: Please report and fix Config.txt manually according to documentation")
                sys.exit(1)
        c+=1
        if c<len(steps) and step != "GENERAL":
            file.write(", ")