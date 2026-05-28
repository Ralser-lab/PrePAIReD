# created on 15.01.26 by Dominik Bierbaum
# error check for correct input variables, contains folling functions:
# read_dataframe:  Check if filepath exists and read the dataframe header
# check_for_existence: Check if value exists in list
# check_for_range: Check if a value is in the range of allowed values
# check_for_path_existence: Check if filepath exists
# check_for_length: Check if l is a list and return its length
# sheet_existence: Find the metadata sheet in an Excel file
# allowed_values_conversion: Convert allowed values string for comparison to user input
# operation_selection: Run the correct operation for every value
# error_check: This function checks parameters and values from the user input for right datatype, default values and allowed ranges for the called modules

import os
import numpy as np
import pandas as pd
from pyarrow.parquet import ParquetFile
import pyarrow as pa

# test data
# metadata_file = "Test_data/Metadata.xlsx"
# data_file = "Test_data/data_00.tsv"

def update_log(logfile, message):
    log = open(logfile, "a")
    log.write(message+"\n")
    log.close()

# support functions
def read_dataframe(filepath: str, log: str, header_only = False, column=None):
    """
    Check if filepath exists and read the dataframe header
    :param filepath: str, path to dataframe
    :param log: str, logfile
    :param header_only: bool, if True only the header will be read
    :column: str or None, column to be read-in
    :return: column names of the dataframe
    """
    # print(filepath)
    
    if header_only:
        nrows = 0
    else:
        nrows = None
        
    if column is not None:
        usecols = [column]
        columns = [column]
    else:
        usecols = None
        columns = None
    
    if not isinstance(filepath, str):  # check if filepath is str
        update_log(log,f"###ERROR###: {filepath} is not a string.") 
        raise TypeError(f"{filepath} is not a string.")
    if filepath.endswith(".xlsx"):  # read excel
        #print(filepath)
        df = pd.read_excel(filepath)  # read only the first 10 lines to save time
    elif filepath.endswith(".txt") or filepath.endswith(".tsv"):  # read tsv files
        df = pd.read_csv(filepath, sep="\t", nrows=nrows, usecols=usecols)  # read only the first line to save time
    elif filepath.endswith(".csv"):  # read csv files
        df = pd.read_csv(filepath, sep=",", nrows=nrows, usecols=usecols)  # read only the first line to save time
    elif filepath.endswith(".parquet"):  # read parquet files
        if columns is not None:
            df = pd.read_parquet(filepath, columns=columns)  
        elif nrows is not None:
            pf = ParquetFile(filepath) 
            first_ten_rows = next(pf.iter_batches(batch_size = 1)) 
            df = pa.Table.from_batches([first_ten_rows]).to_pandas()
        else:
            df = pd.read_parquet(filepath)  
    multiple_columns = df.columns[df.columns.duplicated()].unique()  # check for uniqueness
    if len(multiple_columns) > 0:
        update_log(log,f"###ERROR###: Dataframe columns are not unique in metadata, {list(multiple_columns)} exist(s) multiple "
                         f"times.") 
        raise ValueError(f"Dataframe columns are not unique in metadata, {list(multiple_columns)} exist(s) multiple "
                         f"times.")
    if df.columns.isna().any():  # check for NaN column names
        update_log(log,f"###ERROR###: Dataframe columns contains NaN column names in metadata at position "
                         f"{df.columns[df.columns.isna()]}.") 
        raise ValueError(f"Dataframe columns contains NaN column names in metadata at position "
                         f"{df.columns[df.columns.isna()]}.")
    else:
        return df


def check_for_existence(value, l: list, module: str, variable: str, log: str):
    """
    Check if value exists in list
    :param value: str, int or float
    :param l: list, contains the allowed values
    :param module: str, contains the current module
    :param variable: str, contains the current variable
    :param log: str, logfile
    :return: Nothing if value is in list
    """
    if not isinstance(l, list):  # check if list is list
        update_log(log,f"###ERROR###: {l} is not a list.") 
        raise TypeError(f"{l} is not a list.")
    if str(value) in l:
        pass
    else:
        update_log(log,f"###ERROR###: for variable {variable}: {value} is not in allowed list.") 
        raise ValueError(f"Error in module: {module} for variable: {variable}: {value} is not in allowed list.")


def check_for_range(value, lower_limit, upper_limit, lower_inclusive, upper_inclusive, module: str, variable: str, log: str):
    """
    Check if a value is in the range of allowed values
    :param value: float, needs to be between lower_inclusive and upper_inclusive
    :param lower_limit: int or float, smallest possible value
    :param upper_limit: int or float, biggest possible value
    :param lower_inclusive: str, will be ( (not included) or [ (included)
    :param upper_inclusive: str, will be ) (not included) or ] (included)
    :param module: str, contains the current module
    :param variable: str, contains the current variable
    :param log: str, logfile
    :return: Nothing if all is fine
    """
    try:  # check if value is not a number
        float(value)
        pass
    except:
        update_log(log,f"###ERROR###: for variable {variable}: {value} is not a number.") 
        raise ValueError(f"{value} is not a number.")
    if not isinstance(lower_inclusive, str):  # check if lower_limit is str
        update_log(log,f"###ERROR###: {lower_inclusive} is not a string.") 
        raise TypeError(f"{lower_inclusive} is not a string.")
    if not isinstance(upper_inclusive, str):  # check if upper_limit is str
        update_log(log,f"###ERROR###: {upper_inclusive} is not a string.") 
        raise TypeError(f"{upper_inclusive} is not a string.")

    if lower_inclusive == "(":  # check if value is out of lower range
        if value > lower_limit:
            pass
        else:
            update_log(log,f"###ERROR###: for variable {variable}: {value} is not in allowed range, must be > "
            f"{lower_limit}.") 
            raise ValueError(
                f"Error in module: {module} for variable: {variable}: {value} is not in allowed range, must be > "
                f"{lower_limit}.")
    if lower_inclusive == "[":  # check if value is out of lower range
        if value >= lower_limit:
            pass
        else:
            update_log(log,f"###ERROR###: for variable {variable}: {value} is not in allowed range, must be >= "
            f"{lower_limit}.") 
            raise ValueError(
                f"Error in module: {module} for variable: {variable}: {value} is not in allowed range, must be >= "
                f"{lower_limit}.")
    if upper_inclusive == ")":  # check if value is out of upper range
        if value < upper_limit:
            pass
        else:
            update_log(log,f"###ERROR###: for variable {variable}: {value} is not in allowed range, must be < "
            f"{upper_limit}.") 
            raise ValueError(
                f"Error in module: {module} for variable: {variable}: {value} is not in allowed range, must be < "
                f"{upper_limit}.")
    if upper_inclusive == "]":  # check if value is out of upper range
        if value <= upper_limit:
            pass
        else:
            update_log(log,f"###ERROR###: for variable {variable}: {value} is not in allowed range, "
                             f"must be <= {upper_limit}.") 
            raise ValueError(f"Error in module: {module} for variable: {variable}: {value} is not in allowed range, "
                             f"must be <= {upper_limit}.")


def check_for_path_existence(filepath: str, log: str):
    """
    Check if filepath exists
    :param filepath: str, path to dataframe
    :param log: str, logfile
    :return: Nothing if filepath exists
    """
    if not isinstance(filepath, str):  # check if filepath is str
        update_log(log,f"###ERROR###: File path {filepath} is not a string.") 
        raise TypeError(f"{filepath} is not a string.")
    if os.path.exists(filepath):
        pass
    else:
        update_log(log,f"###ERROR###: File path {filepath} does not exist.") 
        raise ValueError(f"{filepath} does not exist.")


def check_for_length(l: list, module: str, variable: str, log: str):
    """
    Check if l is a list and return its length
    :param l: list, contains the allowed values
    :param module: str, contains the current module
    :param variable: str, contains the current variable
    :param log: str, logfile
    :return: length of l
    """
    if not isinstance(l, list):  # check if list is a list
        update_log(log,f"###ERROR###: for variable {variable}: {l} is not a list.") 
        raise TypeError(f"Error in module: {module} for variable: {variable}: {l} is not a list.")
    return len(l)


def sheet_existence(excel_file: str, sheet_name: str, log: str):
    """
    Find the metadata sheet in an Excel file
    :param excel_file: str, Excel file to check
    :param sheet_name: str, sheet to find in Excel
    :param log: str, logfile
    :return: Nothing if sheet name exists
    """
    if not isinstance(excel_file, str):  # check if excel_file is str
        update_log(log,f"###ERROR###: Excel file {excel_file} is not a string.") 
        raise TypeError(f"{excel_file} is not a string.")
    if not isinstance(sheet_name, str):  # check if sheet_name is str
        update_log(log,f"###ERROR###: Sheet {sheet_name} is not a string.") 
        raise TypeError(f"{sheet_name} is not a string.")
    try:
        xl = pd.read_excel(excel_file, sheet_name=sheet_name)  # read excel
    except:
        update_log(log,f"###ERROR###: Sheet {sheet_name} not found in excel file.") 
        raise ValueError(f"{sheet_name} not found in excel file.")


def allowed_values_conversion(range_options: list, log: str):
    """
    Convert allowed values string for comparison to user input
    :param range_options: allowed range including the boundaries and values
    :param log: str, logfile
    :return: lower_inclusive, upper_inclusive, lower_limit, upper_limit
    """
    allowed_values = range_options[0]
    lower_inclusive = allowed_values[0]  # ( or [
    upper_inclusive = allowed_values[-1]  # ) or ]
    allowed_values = allowed_values[1:-1]  # allowed values without parenthesis
    allowed_values = allowed_values.split(",")
    lower_limit = -np.inf if allowed_values[0].strip() == "-np.inf" else float(allowed_values[0].strip())
    upper_limit = np.inf if allowed_values[1].strip() == "np.inf" else float(allowed_values[1].strip())
    return lower_inclusive, upper_inclusive, lower_limit, upper_limit


def operation_selection(value, target_dict: dict, operation_dictionary: dict, user_dict: dict, operation: str,
                        module: str, variable: str, metadata_file: str, data_file: str, log: str):
    """
    Run the correct operation for every value
    :param value: user input
    :param target_dict: dict, nested dict entries for variable, will be used to identify allowed values
    :param operation_dictionary: dict, operation dict to
    :param user_dict: dict, dictionary with current accepted user parameters
    :param operation: operation to perform
    :param module: str, contains the current module
    :param variable: str, contains the current variable
    :param metadata_file: str, metadata filename
    :param data_file: str, proteomics filename, f.e. DIA-NN output
    :param log: str, logfile
    """
    #print(value)
    current_operation = operation
    if value is not None:  # skip operations if a value is nan
        if current_operation == "sheet existence if xls(x)":
            if metadata_file[:-5] == ".xlsx":
                sheet_existence(metadata_file, value, log=log)
        if current_operation == "existence in metadata.columns":
            metadata_df = read_dataframe(metadata_file, log=log)
            check_for_existence(value, l=list(metadata_df.columns), module=module, variable=variable, log=log)
            if variable in ["metadata_filename_column","metadata_sample_column"]:  # add check for uniqueness and nans in file.name column
                multiple_files = metadata_df.loc[metadata_df[value].duplicated(), value].unique()  # check for uniqueness
                if len(multiple_files) > 0:
                    update_log(log, f"###ERROR###: File or sample names are not unique in metadata, {list(multiple_files)} exist(s) multiple "
                     f"times.") 
                    raise ValueError(
                        f"File or samples names are not unique in metadata, {list(multiple_files)} exist(s) multiple "
                        f"times.")
                if metadata_df[value].isna().any():  # check for NaN column names
                    update_log(log, "###ERROR###: File or sample names contain NaN column names in metadata.") 
                    raise ValueError("File or sample names contain NaN column names in metadata.")
        if current_operation == "existence in metadata.columns and value in metadata[column]":
            if value and len(value) >= 2:
                metadata_df = read_dataframe(metadata_file, log=log)
                check_for_existence(value[0], l=list(metadata_df.columns), module=module, variable=variable, log=log)
                # print(df_cols)
                # print(value[0])
                series = metadata_df[value[0]].tolist()
                #series = series[value[0]].tolist()
                # print(series)
                check_for_existence(value[1], l=series, module=module, variable=variable, log=log)
        if current_operation == "existence in metadata_df.columns (; separated)":
            metadata_df = read_dataframe(metadata_file, log=log)
            values = value.split(";")
            for val in values:
                check_for_existence(val, l=list(metadata_df.columns), module=module, variable=variable, log=log)
        if current_operation == "elements existence in metadata.columns":
            metadata_df = read_dataframe(metadata_file, log=log)
            for val in value:
                check_for_existence(val, l=list(metadata_df.columns), module=module, variable=variable, log=log)
        if current_operation == "value in metadata[sample_type_column]":
            metadata_df = read_dataframe(metadata_file, log=log)
            check_for_existence(value, l=list(metadata_df[user_dict["sample_type_column"]]), module=module, variable=variable, log=log)
        if current_operation == "value in metadata[combat_fitsamples_meta_column]":
            metadata_df = read_dataframe(metadata_file, log=log)
            check_for_existence(value, l=list(metadata_df[user_dict["combat_fitsamples_meta_column"]]), module=module, variable=variable, log=log)
        if current_operation == "value in metadata[adjust_fitsamples_meta_column]":
            metadata_df = read_dataframe(metadata_file, log=log)
            check_for_existence(value, l=list(metadata_df[user_dict["adjust_fitsamples_meta_column"]]), module=module, variable=variable, log=log)
        if current_operation == "existence in data_df.columns":
            if variable == "input_filename_column":  # add check for uniqueness and nans in file.name column
                # Check wheter column exists
                data_df = read_dataframe(data_file, log=log, header_only=True)
                check_for_existence(value, l=list(data_df.columns), module=module, variable=variable, log=log)
                #multiple_files = data_df.loc[data_df[value].duplicated(), value].unique()  # check for uniqueness
                #if len(multiple_files) > 0:
                #    update_log(log, f"###ERROR###: File names are not unique in data, {list(multiple_files)} exist(s) multiple "
                #    f"times.") 
                #    raise ValueError(
                #        f"File names are not unique in data, {list(multiple_files)} exist(s) multiple "
                #        f"times.")
                # read in only respective column
                data_df = read_dataframe(data_file, log=log, column=value)
                if data_df[value].isna().any():  # check for NaN column names
                    raise ValueError("File names contain NaN column names in data.")
                # if file names in metadata and data are unique and without nans, check that file names are same in
                # both files
                metadata_df = read_dataframe(metadata_file, log=log)
                metadata_set = set(metadata_df[user_dict["metadata_filename_column"]])
                data_set = set(data_df[user_dict["input_filename_column"]])
                #only_in_metadata = metadata_set - data_set
                only_in_data = data_set - metadata_set
                if only_in_data:
                    #update_log(log, f"###ERROR###: Differences found:\n"
                    #                 #f"Only in metadata: {list(only_in_metadata)}\n"
                    #                 f"Only in data: {list(only_in_data)}")
                    raise ValueError(f"Differences found:\n"
                                     #f"Only in metadata: {list(only_in_metadata)}\n"
                                     f"Only in data: {list(only_in_data)}")
            else:
                data_df = read_dataframe(data_file, log=log, header_only=True)
                check_for_existence(value, l=list(data_df.columns), module=module, variable=variable, log=log)

        if current_operation == "existence in data_df.columns (; separated)":
            data_df = read_dataframe(data_file, log=log, header_only=True)
            values = value.split(";")
            for val in values:
                check_for_existence(val, l=list(data_df.columns), module=module, variable=variable, log=log)
        if current_operation == "existence in input_ids_columns":
            # print(user_dict)
            input_ids_columns = user_dict["input_ids_columns"]
            #print(input_ids_columns)
            input_ids_columns = input_ids_columns.split(";")
            #print(input_ids_columns)
            check_for_existence(value, l=input_ids_columns, module=module, variable=variable, log=log)
        if current_operation == "existence in input_ids_columns or Proteotypic":
            # print(user_dict)
            input_ids_columns = user_dict["input_ids_columns"]
            #print(input_ids_columns)
            input_ids_columns = input_ids_columns.split(";")
            input_ids_columns.append("Proteotypic")
            input_ids_columns.append("include_correction")
            check_for_existence(value, l=input_ids_columns, module=module, variable=variable, log=log)
        if current_operation == "existence of path":
            check_for_path_existence(filepath=value, log=log)
        if current_operation == "check":
            check_for_existence(value, l=target_dict["Options"], module=module, variable=variable, log=log)
        if current_operation == "check len(adjust_cols) and el in AllowedValues":
            adjust_cols = user_dict["adjust_cols"]  # get the list with column names
            allowed_values = operation_dictionary["correct_method"]["adjust_confounders"]["adjust_cols_dtypes"][
                "Options"]  # dtypes
            length_cols = check_for_length(adjust_cols, module=module, variable=variable, log=log)
            length_dtypes = check_for_length(value, module=module, variable=variable, log=log)
            if length_cols != length_dtypes:
                update_log(log, f"###ERROR###: for variable {variable}: Length of dtypes ({length_dtypes}) "
                                 f"and adjusted columns ({length_cols}) is different")
                raise ValueError(f"Error in module: {module} for variable: {variable}: Length of dtypes ({length_dtypes}) "
                                 f"and adjusted columns ({length_cols}) is different")
            for val in value:
                check_for_existence(val, l=allowed_values, module=module, variable=variable, log=log)   
        if current_operation == "check len(combat_cols_adjust_meta) and el in AllowedValues if existence of combat_cols_adjust_meta":   
            if user_dict["combat_cols_adjust_meta"]:
                adjust_cols = user_dict["combat_cols_adjust_meta"]  # get the list with column names
                allowed_values = operation_dictionary["correct_method"]["combat"]["combat_cols_adjust_dtypes"][
                    "Options"]  # dtypes
                length_cols = check_for_length(adjust_cols, module=module, variable=variable, log=log)
                length_dtypes = check_for_length(value, module=module, variable=variable, log=log)
                if length_cols != length_dtypes:
                    update_log(log, f"###ERROR###: for variable {variable}: Length of dtypes ({length_dtypes}) "
                                     f"and adjusted columns ({length_cols}) is different")
                    raise ValueError(f"Error in module: {module} for variable: {variable}: Length of dtypes ({length_dtypes}) "
                                     f"and adjusted columns ({length_cols}) is different")
                for val in value:
                    check_for_existence(val, l=allowed_values, module=module, variable=variable, log=log)              
        if current_operation == "check elements":
            values = value.split(";")
            for val in values:
                check_for_existence(val, l=target_dict["Options"], module=module, variable=variable, log=log)
        if current_operation == "check range":
            allowed_values = target_dict["Options"]
            lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values, log=log)
            check_for_range(value=value, lower_limit=lower_limit, upper_limit=upper_limit,
                            lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                            variable=variable, log=log)
        if current_operation == "check len 2 and range":
            allowed_values = target_dict["Options"]
            allowed_values = allowed_values[0].split("in ")
            allowed_values = allowed_values[-1]
            allowed_values = [allowed_values.strip()]
            lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values, log=log)
            if check_for_length(value, module=module, variable=variable, log=log) != 2:
                update_log(log, f"###ERROR###: for variable {variable}: {value} is supposed to have only 2 elements!")
                raise ValueError(
                    f"Error in module: {module} for variable: {variable}: {value} is supposed to have only 2 elements!")
            for val in value:
                val = float(val)
                check_for_range(value=val, lower_limit=lower_limit, upper_limit=upper_limit,
                                lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                                variable=variable, log=log)
        if current_operation == "check len 3 and range":
            allowed_values = target_dict["Options"]
            allowed_values = allowed_values[0].split("in ")
            allowed_values = allowed_values[-1]
            allowed_values = [allowed_values.strip()]
            lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values, log=log)
            if check_for_length(value, module=module, variable=variable, log=log) != 3:
                update_log(log, f"###ERROR###: for variable {variable}: {value} is supposed to have only 3 elements!")
                raise ValueError(
                    f"Error in module: {module} for variable: {variable}: {value} is supposed to have only 3 elements!")
            for val in value:
                val = float(val)
                check_for_range(value=val, lower_limit=lower_limit, upper_limit=upper_limit,
                                lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                                variable=variable, log=log)
        if current_operation == "check elements (; or _ separated)":
            values = str(value).split(";")
            values = [el for item in values for el in item.split("_")]
            allowed_values = target_dict["Options"]
            lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values, log=log)
            for val in values:
                val = float(val)
                check_for_range(value=val, lower_limit=lower_limit, upper_limit=upper_limit,
                                lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                                variable=variable, log=log)
        if current_operation == "check range element-wise":
            allowed_values = target_dict["Options"]
            allowed_values = allowed_values[0].split("in ")
            allowed_values = allowed_values[-1]
            allowed_values = [allowed_values.strip()]
            lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values, log=log)
            for val in value:
                val = float(val)
                check_for_range(value=val, lower_limit=lower_limit, upper_limit=upper_limit,
                                lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                                variable=variable, log=log)

    else:  # check if NoneType is allowed
        allowed_types = target_dict["Type"]
        allowed_types = allowed_types.split(";")
        if pd.isna(target_dict["Default"]):
            pass
        else:
            if "None" not in allowed_types:
                update_log(log, f"###ERROR###: for variable {variable}: NoneType is not allowed.")
                raise TypeError(f"Error in module {module} for variable {variable}: NoneType is not allowed.")

# main function #
# read DIANN file header, metadata
def error_check(value, variable: str, module: str, flag_dictionary: dict, operation_dictionary: dict, user_dict: dict,
                metadata_file: str, data_file: str, log: str, parameters_list: list,
                current_flag=None, current_flag_value=None):
    """
    This function checks parameters and values from the user input for right datatype, default
    values and allowed ranges for the called modules.

    :param value: any dtype, input from user to check for any errors
    :param variable: str, variable to check value for
    :param module: str, module name
    :param flag_dictionary: dictionary containing flag informations
    :param operation_dictionary: if a variable is associated to a flag, identify here via flag and flag_value
    :param current_flag: str, flag is needed to identify associated parameters
    :param current_flag_value: str, flag_value is needed to identify associated parameters
    :param user_dict: dictionary with current accepted user parameters
    :param metadata_file: str, metadata filename
    :param data_file: str, proteomics filename, f.e. DIA-NN output
    :param log: str, logfile
    :param parameters_list: list, contains all user-derived parameters to cross-check with mandatory values
    :return current_flag, current_flag_value, str, both is needed to identify certain variables in operation dict
    """

    # check all arguments for right type
    # if not isinstance(variable, str):  # check if variable is really a str
    #    raise TypeError(f"{variable} is not a string.")
    # if not isinstance(module, str):  # check if module is really a str
    #    raise TypeError(f"{module} is not a string.")
    # if not isinstance(flag_dictionary, dict):  # check if dict is really a dict
    #    raise TypeError(f"{flag_dictionary} is not a dictionary.")
    # if not isinstance(operation_dictionary, dict):  # check if dict is really a dict
    #    raise TypeError(f"{operation_dictionary} is not a dictionary.")

    # identify value via module and parameter and compare to dict hit
    # make a try except statement:
    # try: store flag as flag_var and flag_value as flag_value_var for temporary variables, example input_origin = Flag,
    # diann = Flag_value, get operations from flag_dict
    # except: if no flag and flag_value is found, use temporary flag and value to check in parameter dict

    if variable in flag_dictionary[module].keys():
        # print(f"{variable} is in Flag dictionary")
        # print(f"Latest flag: {current_flag}, latest flag value: {current_flag_value}, "
        #       f"variable: {variable}, value: {value}, module: {module}")
        # print(f"current parameters list: {parameters_list}")
        if current_flag is None:  # exception for first run
            pass
        else:
            if len(parameters_list) == 0:  # look only at parameters which are present
                pass
            else:
                if len([key for key in operation_dictionary[current_flag].keys()]) > 1:  # some flag values don't need to be checked before there is a parameter
                    if ";" in current_flag_value:
                        flag_values = current_flag_value.split(";")
                        for flag_value in flag_values:
                            mandatory_parameters = [key for key in operation_dictionary[current_flag][flag_value].keys()]  # get mandatory parameters
                            mandatory_parameters.remove("Language")  # Language is not used in analysis
                            # print(f" mandatory parameters: {mandatory_parameters}")
                            if set(mandatory_parameters).issubset(set(parameters_list)):  # check if mandatory parameters are in list
                                # print("All parameters found.")
                                # parameters_list = []
                            else:
                                missing = set(mandatory_parameters) - set(parameters_list)  # get missing prarameters
                                # print(f"Parameters missing: {missing}")
                                update_log(log,
                                           f"###ERROR###: Missing parameters for {current_flag} and {current_flag_value}: {missing}")
                                raise ValueError(f"###ERROR###: Missing parameters for {current_flag} and {current_flag_value}: "
                                                 f"{missing}")
                    else:
                        mandatory_parameters = [key for key in operation_dictionary[current_flag][current_flag_value].keys()]
                        mandatory_parameters.remove("Language")
                        # print(f" mandatory parameters: {mandatory_parameters}")
                        if set(mandatory_parameters).issubset(set(parameters_list)):
                            # print("All parameters found.")
                            # parameters_list = []
                        else:
                            missing = set(mandatory_parameters) - set(parameters_list)
                            # print(f"Parameters missing: {missing}")
                            update_log(log,
                                       f"###ERROR###: Missing parameters for {current_flag} and {current_flag_value}: {missing}")
                            raise ValueError(
                                f"###ERROR###: Missing parameters for {current_flag} and {current_flag_value}: {missing}")
        dict_to_use = flag_dictionary[module][variable]  # get the current flag information
        # print(dict_to_use)
        current_flag = variable
        current_flag_value = str(value)

    else:
        # print(f"{variable} is in Operation dictionary")
        if current_flag_value is None:
            update_log(log, "###ERROR###: Flag values can not be None.")
            raise ValueError("Flag values can not be None.")
        if ";" in current_flag_value:  # make exceptions if you have to loop over several flag values seperated by ;
            flag_values = current_flag_value.split(";")
            # dict_to_use = operation_dictionary[current_flag][flag_values[0]][variable]
            for flag_value in flag_values:
                if variable in operation_dictionary[current_flag][flag_value]:
                    # print(flag_value, variable)
                    dict_to_use = operation_dictionary[current_flag][flag_value][variable]
                    parameters_list.append(variable)
                    # print(value, type(value))
                    # print(dict_to_use["Default"], type(dict_to_use["Default"]))

                    current_operation = dict_to_use["Check Operation"]  # get the current operation to do
                    if str(value).strip() == dict_to_use[
                        "Default"] and "existence" not in current_operation:  # check for default
                        # print(f"{value} is default value.")
                        pass
                    else:
                        # print(f"{value} is not default")
                        # current_operation = dict_to_use["Check Operation"]  # get the current operation to do
                        # print(current_operation)
                        operation_selection(value=value, target_dict=dict_to_use,
                                            operation_dictionary=operation_dictionary,
                                            operation=current_operation, module=module, variable=variable,
                                            user_dict=user_dict, metadata_file=metadata_file, data_file=data_file,
                                            log=log)
                else:
                    pass
        else:
            try:
                dict_to_use = operation_dictionary[current_flag][current_flag_value][variable]
                parameters_list.append(variable)
            except:
                update_log(log,
                           f"###ERROR###: {variable} is not present for flag {current_flag}={current_flag_value}. Please remove manually to avoid any unwanted inconsistencies later on.")
                # print(dict_to_use)

    # print(value, type(value))
    # print(dict_to_use["Default"], type(dict_to_use["Default"]))

    current_operation = dict_to_use["Check Operation"]
    if str(value) == str(dict_to_use["Default"]) and "existence" not in current_operation:  # check for default
        if str(value) == "None" and not "None" in dict_to_use["Type"].split(";"):
            update_log(log, f"###ERROR###: {variable} can not be None.")
            raise ValueError(f"{variable} can not be None.")
        else:
            # print(f"{value} is default value.")
            pass

    else:
        if str(value) == "None" and not "None" in dict_to_use["Type"].split(";"):
            update_log(log, f"###ERROR###: {variable} can not be None.")
            raise ValueError(f"{variable} can not be None.")
        # print(f"{value} is not default")
        # current_operation = dict_to_use["Check Operation"]  # get the current operation to do
        # print(current_operation)
        operation_selection(value=value, target_dict=dict_to_use, operation_dictionary=operation_dictionary,
                            operation=current_operation, module=module, variable=variable, user_dict=user_dict,
                            metadata_file=metadata_file, data_file=data_file, log=log)

    # print(current_flag, current_flag_value)
    # print("Works so far")

    return current_flag, current_flag_value, parameters_list

