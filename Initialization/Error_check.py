# created on 15.01.26 by Dominik Bierbaum
# error check for correct input variables

import os
import numpy as np
import pandas as pd


# test data
# metadata_file = "Test_data/Metadata.xlsx"
# data_file = "Test_data/data_00.tsv"

# support functions
def read_dataframe_header(filepath: str):
    """
    Check if filepath exists and read the dataframe header
    :param filepath: str, path to dataframe
    :return: column names of the dataframe
    """
    if not isinstance(filepath, str):  # check if filepath is str
        raise TypeError(f"{filepath} is not a string.")
    if filepath[:-5] == ".xlsx":  # read excel
        df = pd.read_excel(filepath)  # read only the first 10 lines to save time
    elif filepath[:-4] == ".txt" or filepath[:-4] == ".tsv":  # read tsv files
        df = pd.read_csv(filepath, sep="\t")  # read only the first 10 lines to save time
    elif filepath[-4:] == ".csv":  # read csv files
        df = pd.read_csv(filepath, sep=",")  # read only the first 10 lines to save time

    multiple_columns = df.columns[df.columns.duplicated()].unique()  # check for uniqueness
    if len(multiple_columns) > 0:
        raise ValueError(f"Dataframe columns are not unique in metadata, {list(multiple_columns)} exist(s) multiple "
                         f"times.")
    if df.columns.isna().any():  # check for NaN column names
        raise ValueError(f"Dataframe columns contains NaN column names in metadata at position "
                         f"{df.columns[df.columns.isna()]}.")
    else:
        return df


def check_for_existence(value, l: list, module: str, variable: str):
    """
    Check if value exists in list
    :param value: str, int or float
    :param l: list, contains the allowed values
    :param module: str, contains the current module
    :param variable: str, contains the current variable

    :return: Nothing if value is in list
    """
    if not isinstance(l, list):  # check if list is list
        raise TypeError(f"{l} is not a list.")
    if value in l:
        pass
    else:
        raise ValueError(f"Error in module: {module} for variable: {variable}: {value} is not in allowed list.")


def check_for_range(value, lower_limit, upper_limit, lower_inclusive, upper_inclusive, module: str, variable: str):
    """
    check if a value is in the range of allowed values
    :param value: float, needs to be between lower_inclusive and upper_inclusive
    :param lower_limit: int or float, smallest possible value
    :param upper_limit: int or float, biggest possible value
    :param lower_inclusive: str, will be ( (not included) or [ (included)
    :param upper_inclusive: str, will be ) (not included) or ] (included)
    :param module: str, contains the current module
    :param variable: str, contains the current variable

    :return: Nothing if all is fine
    """
    try:  # check if value is not a number
        float(value)
        pass
    except:
        raise ValueError(f"{value} is not a number.")
    if not isinstance(lower_inclusive, str):  # check if lower_limit is str
        raise TypeError(f"{lower_inclusive} is not a string.")
    if not isinstance(upper_inclusive, str):  # check if upper_limit is str
        raise TypeError(f"{upper_inclusive} is not a string.")

    if lower_inclusive == "(":  # check if value is out of lower range
        if value > lower_limit:
            pass
        else:
            raise ValueError(
                f"Error in module: {module} for variable: {variable}: {value} is not in allowed range, must be > "
                f"{lower_limit}.")
    if lower_inclusive == "[":  # check if value is out of lower range
        if value >= lower_limit:
            pass
        else:
            raise ValueError(
                f"Error in module: {module} for variable: {variable}: {value} is not in allowed range, must be >= "
                f"{lower_limit}.")
    if upper_inclusive == ")":  # check if value is out of upper range
        if value < upper_limit:
            pass
        else:
            raise ValueError(
                f"Error in module: {module} for variable: {variable}: {value} is not in allowed range, must be < "
                f"{upper_limit}.")
    if upper_inclusive == "]":  # check if value is out of upper range
        if value <= upper_limit:
            pass
        else:
            raise ValueError(f"Error in module: {module} for variable: {variable}: {value} is not in allowed range, "
                             f"must be <= {upper_limit}.")


def check_for_path_existence(filepath: str):
    """
    Check if filepath exists
    :param filepath: str, path to dataframe
    :return: Nothing if filepath exists
    """
    if not isinstance(filepath, str):  # check if filepath is str
        raise TypeError(f"{filepath} is not a string.")
    if os.path.exists(filepath):
        pass
    else:
        raise ValueError(f"{filepath} does not exist.")


def check_for_length(l: list, module: str, variable: str):
    """
    Check if l is a list and return its length
    :param l: list, contains the allowed values
    :param module: str, contains the current module
    :param variable: str, contains the current variable
    :return: length of l
    """
    if not isinstance(l, list):  # check if list is a list
        raise TypeError(f"Error in module: {module} for variable: {variable}: {l} is not a list.")
    return len(l)


def sheet_existence(excel_file: str, sheet_name: str):
    """
    Find the metadata sheet in an Excel file
    :param excel_file: str, Excel file to check
    :param sheet_name: str, sheet to find in Excel
    :return: Nothing if sheet name exists
    """
    if not isinstance(excel_file, str):  # check if excel_file is str
        raise TypeError(f"{excel_file} is not a string.")
    if not isinstance(sheet_name, str):  # check if sheet_name is str
        raise TypeError(f"{sheet_name} is not a string.")
    try:
        xl = pd.read_excel(excel_file, sheet_name=sheet_name)  # read excel
    except:
        raise ValueError(f"{sheet_name} not found in excel file.")


def allowed_values_conversion(range_options: list):
    """
    Convert allowed values string for comparison to user input
    :param range_options: allowed range including the boundaries and values
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
                        module: str, variable: str, metadata_file: str, data_file: str):
    """
    run the correct operation for every value
    :param value: user input
    :param target_dict: dict, nested dict entries for variable, will be used to identify allowed values
    :param operation_dictionary: dict, operation dict to
    :param user_dict: dict, dictionary with current accepted user parameters
    :param operation: operation to perform
    :param module: str, contains the current module
    :param variable: str, contains the current variable
    :param metadata_file: str, metadata filename
    :param data_file: str, proteomics filename, f.e. DIA-NN output
    """
    print(operation)
    print(value)
    current_operation = operation

    if current_operation == "sheet existence if xls(x)":
        if metadata_file[:-5] == ".xlsx":
            sheet_existence(metadata_file, value)
    if current_operation == "existence in metadata.columns":
        metadata_df = read_dataframe_header(metadata_file)
        check_for_existence(value, l=list(metadata_df.columns), module=module, variable=variable)
        if variable == "metadata_filename_column":  # add check for uniqueness and nans in file.name column
            multiple_files = metadata_df.loc[df[value].duplicated(), value].unique()  # check for uniqueness
            if len(multiple_files) > 0:
                raise ValueError(
                    f"File names are not unique in metadata, {list(multiple_files)} exist(s) multiple "
                    f"times.")
            if metadata_df[value].isna().any():  # check for NaN column names
                raise ValueError(f"File names contain NaN column names in metadata.")
    if current_operation == "existence in metadata.columns and value in metadata[column]":
        metadata_df = read_dataframe_header(metadata_file)
        check_for_existence(value[0], l=list(metadata_df.columns), module=module, variable=variable)
        # print(df_cols)
        # print(value[0])
        series = pd.read_excel(metadata_file, usecols=[value[0]])
        series = series[value[0]].tolist()
        # print(series)
        check_for_existence(value[1], l=series, module=module, variable=variable)
    if current_operation == "existence in metadata_df.columns (; separated)":
        metadata_df = read_dataframe_header(metadata_file)
        values = value.split(";")
        for val in values:
            check_for_existence(val, l=list(metadata_df.columns), module=module, variable=variable)
    if current_operation == "elements existence in metadata.columns":
        metadata_df = read_dataframe_header(metadata_file)
        for val in value:
            check_for_existence(val, l=list(metadata_df.columns), module=module, variable=variable)
    if current_operation == "existence in data_df.columns":
        data_df = read_dataframe_header(data_file)
        check_for_existence(value, l=list(data_df.columns), module=module, variable=variable)
        if variable == "input_filename_column":  # add check for uniqueness and nans in file.name column
            multiple_files = data_df.loc[df[value].duplicated(), value].unique()  # check for uniqueness
            if len(multiple_files) > 0:
                raise ValueError(
                    f"File names are not unique in data, {list(multiple_files)} exist(s) multiple "
                    f"times.")
            if data_df[value].isna().any():  # check for NaN column names
                raise ValueError(f"File names contain NaN column names in data.")
            # if file names in metadata and data are unique and without nans, check that file names are same in
            # both files
            metadata_set = set(metadata_df[user_dict["metadata_filename_column"]])
            data_set = set(data_df[user_dict["input_filename_column"]])
            only_in_metadata = metadata_set - data_set
            only_in_data = data_set - metadata_set
            if only_in_metadata or only_in_data:
                raise ValueError(f"Differences found:\n"
                                 f"Only in metadata: {list(only_in_metadata)}\n"
                                 f"Only in data: {list(only_in_data)}")

    if current_operation == "existence in data_df.columns (; separated)":
        data_df = read_dataframe_header(data_file)
        values = value.split(";")
        for val in values:
            check_for_existence(val, l=list(data_df.columns), module=module, variable=variable)
    if current_operation == "existence in input_ids_columns":
        # print(user_dict)
        input_ids_columns = user_dict["input_ids_columns"]
        print(input_ids_columns)
        input_ids_columns = input_ids_columns.split(";")
        print(input_ids_columns)
        check_for_existence(value, l=input_ids_columns, module=module, variable=variable)
    if current_operation == "existence in input_ids_columns or Proteotypic":
        # print(user_dict)
        input_ids_columns = user_dict["input_ids_columns"]
        print(input_ids_columns)
        input_ids_columns = input_ids_columns.split(";")
        input_ids_columns.append("Proteotypic")
        check_for_existence(value, l=input_ids_columns, module=module, variable=variable)
    if current_operation == "existence of path":
        check_for_path_existence(filepath=value)
    if current_operation == "check":
        check_for_existence(value, l=target_dict["Options"], module=module, variable=variable)
    if current_operation == "check len(adjust_cols) and el in AllowedValues":
        adjust_cols = user_dict["adjust_cols"]  # get the list with column names
        allowed_values = operation_dictionary["correct_method"]["adjust_confounders"]["adjust_cols_dtypes"][
            "Options"]  # dtypes
        length_cols = check_for_length(adjust_cols, module=module, variable=variable)
        length_dtypes = check_for_length(value, module=module, variable=variable)
        if length_cols != length_dtypes:
            raise ValueError(f"Error in module: {module} for variable: {variable}: Length of dtypes ({length_dtypes}) "
                             f"and adjusted columns ({length_cols}) is different")
        for val in value:
            check_for_existence(val, l=allowed_values, module=module, variable=variable)
    if current_operation == "check elements":
        values = value.split(";")
        for val in values:
            check_for_existence(val, l=target_dict["Options"], module=module, variable=variable)
    if current_operation == "check range":
        allowed_values = target_dict["Options"]
        lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values)
        check_for_range(value=value, lower_limit=lower_limit, upper_limit=upper_limit,
                        lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                        variable=variable)
    if current_operation == "check len 2 and range":
        allowed_values = target_dict["Options"]
        lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values)
        if check_for_length(value, module=module, variable=variable) != 2:
            raise ValueError(
                f"Error in module: {module} for variable: {variable}: {value} is supposed to have only 2 elements!")
        for val in value:
            val = float(val)
            check_for_range(value=val, lower_limit=lower_limit, upper_limit=upper_limit,
                            lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                            variable=variable)
    if current_operation == "check len 3 and range":
        allowed_values = target_dict["Options"]
        lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values)
        if check_for_length(value, module=module, variable=variable) != 3:
            raise ValueError(
                f"Error in module: {module} for variable: {variable}: {value} is supposed to have only 3 elements!")
        for val in value:
            val = float(val)
            check_for_range(value=val, lower_limit=lower_limit, upper_limit=upper_limit,
                            lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                            variable=variable)
    if current_operation == "check elements (; or _ separated)":
        values = value.split(";").split("_")
        allowed_values = target_dict["Options"]
        lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values)
        for val in values:
            val = float(val)
            check_for_range(value=val, lower_limit=lower_limit, upper_limit=upper_limit,
                            lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                            variable=variable)
    if current_operation == "check range element-wise":
        allowed_values = target_dict["Options"]
        allowed_values = allowed_values[0].split("in ")
        allowed_values = allowed_values[-1]
        allowed_values = [allowed_values.strip()]
        lower_inclusive, upper_inclusive, lower_limit, upper_limit = allowed_values_conversion(allowed_values)
        for val in value:
            val = float(val)
            check_for_range(value=val, lower_limit=lower_limit, upper_limit=upper_limit,
                            lower_inclusive=lower_inclusive, upper_inclusive=upper_inclusive, module=module,
                            variable=variable)


# main function #
# read DIANN file header, metadata
def error_check(value, variable: str, module: str, flag_dictionary: dict, operation_dictionary: dict, current_flag: str,
                current_flag_value: str, user_dict: dict, metadata_file: str, data_file: str):
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
    :return current_flag, current_flag_value, str, both is needed to identify certain variables in operation dict
    """

    # check all arguments for right type
    if not isinstance(variable, str):  # check if variable is really a str
        raise TypeError(f"{variable} is not a string.")
    if not isinstance(module, str):  # check if module is really a str
        raise TypeError(f"{module} is not a string.")
    if not isinstance(flag_dictionary, dict):  # check if dict is really a dict
        raise TypeError(f"{flag_dictionary} is not a dictionary.")
    if not isinstance(operation_dictionary, dict):  # check if dict is really a dict
        raise TypeError(f"{operation_dictionary} is not a dictionary.")

    # identify value via module and parameter and compare to dict hit
    # make a try except statement:
    # try: store flag as flag_var and flag_value as flag_value_var for temporary variables, example input_origin = Flag,
    # diann = Flag_value, get operations from flag_dict
    # except: if no flag and flag_value is found, use temporary flag and value to check in parameter dict
    if variable in flag_dictionary[module].keys():
        print(f"{variable} is in Flag dictionary")
        dict_to_use = flag_dictionary[module][variable]  # get the current flag information
        print(dict_to_use)
        current_flag = variable
        current_flag_value = str(value)

    else:
        print(f"{variable} is in Operation dictionary")
        if ";" in current_flag_value:  # make exceptions if you have to loop over several flag values seperated by ;
            flag_values = current_flag_value.split(";")
            # dict_to_use = operation_dictionary[current_flag][flag_values[0]][variable]
            for flag_value in flag_values:
                if variable in operation_dictionary[current_flag][flag_value]:
                    print(flag_value, variable)
                    dict_to_use = operation_dictionary[current_flag][flag_value][variable]
                    print(value, type(value))
                    print(dict_to_use["Default"], type(dict_to_use["Default"]))

                    if str(value).strip() == dict_to_use["Default"]:  # check for default
                        print(f"{value} is default value.")
                        pass
                    else:
                        print(f"{value} is not default")
                        current_operation = dict_to_use["Check Operation"]  # get the current operation to do
                        # print(current_operation)
                        operation_selection(value=value, target_dict=dict_to_use,
                                            operation_dictionary=operation_dictionary,
                                            operation=current_operation, module=module, variable=variable,
                                            user_dict=user_dict, metadata_file=metadata_file, data_file=data_file)
                else:
                    pass
        else:
            dict_to_use = operation_dictionary[current_flag][current_flag_value][variable]
            print(dict_to_use)

    print(value, type(value))
    print(dict_to_use["Default"], type(dict_to_use["Default"]))

    if str(value) == dict_to_use["Default"]:  # check for default
        print(f"{value} is default value.")
        pass
    elif pd.isna(value and pd.isna(dict_to_use["Default"])):  # check for None entries
        print(f"{value} is default value.")
        pass
    else:
        print(f"{value} is not default")
        current_operation = dict_to_use["Check Operation"]  # get the current operation to do
        # print(current_operation)
        operation_selection(value=value, target_dict=dict_to_use, operation_dictionary=operation_dictionary,
                            operation=current_operation, module=module, variable=variable, user_dict=user_dict,
                            metadata_file=metadata_file, data_file=data_file)

    print(current_flag, current_flag_value)
    print("Works so far")

    return current_flag, current_flag_value


