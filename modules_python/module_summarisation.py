"""
Input: wide-format data
Output: Summarised wide-format data

    - print_args: decorator that prints arguments of a function (sanity checks)
    - calc_ref_intensity: calculates reference intensities
    - apply_protein_filter: filter proteins by their number of runs that they are detected in.
    - remove_rows_cols: remove rows and columns from a matrix.
    - get_peptide_ratios: calculate peptide ratios between runs
    - get_reference_intensities: calculate reference intensities for each protein in each run
    - estimate_protein_profile: estimate protein quantities in each run
    - weighted_harmonic_mean: calculates weighted harmonic mean
    - geometric_mean: calculates the geometric mean of an array x
    - summarize_peptides: summarizes peptide intensities to protein intensities with own approach using a pythagorean mean on representative precursors for a protein
    - transform_profiles_to_longformat_table: reshapes wide format profiles with IDs as rows, samples as columns and abundances as values to long-format with ID, sample and values as columns
    - MaxLFQ_summarization: summarize precursors to protein intensities using MaxLFQ
    - precursor_sum: sums precursor intensities to obtain protein quantities

"""

# TODO
# Adjust doc-strings

#general
import sys
import pandas as pd
import numpy as np
import time
#from pathlib import Path
import os

#own scripts
import module_helper as helper 

#import warnings
#warnings.filterwarnings("ignore")

def weighted_harmonic_mean(x,w=None):
    """
    calculate weighted harmonic mean.
    
    author: BZ
    
    Parameters
    ----------
    x : np.array or list
        Array of floats for which the weighted harmonic mean should be calculated.
    w : np.array or list, optional
        weights associated with values in x for the calculation of the weighted harmonic mean. The default is None, resulting in no weightingfor x.

    Returns
    -------
    float : weighted harmonic mean of x

    """
    x=np.array(x)
    if w is not None:
        w=np.array(w)
    else:
        w=np.ones(len(x))
    return(np.sum(w)/np.sum(w/x))
    


def geometric_mean(x):
    """
    Calculate the geometric mean of an array x.
    
    author: BZ
    
    Parameters
    ----------
    x : np.array
        Array of floats.

    Returns
    -------
    float

    """
    return(np.prod(x)**(1/len(x)))
    


def topN_summarization(data_df, id_mapping_df, filter_mask, module_output, top_n_peptides=3,weights=None,pooling="mean",min_peptides=1, seq_col="Precursor.ID", protein_col = "Protein.Group"):
    """
    summarize peptide intensities to protein intensities
    top_n_peptides most abundant, highly correlated and least imputed precursors over samples are summarized to a protein intensity
    #ranks precursors on number of measured peptides (not imputed ones)
    #ranks precursors on abundance
    #ranks precursors on correlation
    all three ranks are normalized, their harmonic mean is calculated and the average as defined by pooling over top_n_peptides precursors is calculated to retrieve the protein intensity.
    pooling is how to summarize selected precursors to protein group profiles, one of pythagorean means:  {"mean","hmean","gmean"}
    
    author: BZ
    
    Parameters
    ----------
    data_df : pd.DataFrame (old: peptide_profiles)
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    id_mapping_df: pd.DataFrame
        information on mapping from precursor to protein (group) level. rows denote individual precursors
    filter_mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.
    module_output: str
        log string to be filled.
    top_n_peptides : int, optional
        Number of precursors to summarize to protein intensities. The default is 3.
    weights : iterable, optional
        Weights for calculation of harmonic mean over the scores. The default is None, assigning equal weights.
    pooling : str, optional
        one of {"mean","hmean","gmean"} to summarize selected precursors to proteins. "hmean" equals harmonic mean, "gmean" equals geometric mean and "mean" the arithmetic mean. The default is "mean".
    min_peptides : int, optional
        require at least min_peptides precursors to calculate a protein intensity, else the protein intensity is set to 0. The default is 1.
    seq_col : str
        name of index level containing precursor Ids. The default is 'Precursor.Id'.
    protein_col : str
        name of index level containing protein Ids. The default is 'Protein.Group'.

    Returns
    -------
    pg_profiles : pd.DataFrame
        rows equal protein groups, columns correspond to samples and values in the dataframe to intensities.
    id_mapping_df_c: pd.DataFrame
        information on mapping from precursor to protein (group) level. rows denote individual precursors. Index set to seq_col.
    module_output: str
        log string to be filled.

    """
    
    pg_profiles_df=pd.DataFrame(columns=data_df.columns)
    
    if weights is None:
        weights=[1,1,1]

    # unlog data
    data_df.replace(np.nan,0,inplace=True)
        
    module_output += time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tSummarizing precursor quantities to Protein.Group intensities using own method...\n"
    module_output += "\t\t\t\t- summarize top "+str(top_n_peptides)+" precursors by their "+pooling+"\n"
    module_output += "\t\t\t\t- for the summarization precursors were selected by high abundance, correlation and low number of imputed values with following weights: "+", ".join([str(x) for x in weights])+"\n"

    indices=data_df.index
        
    filter_mask_z = filter_mask.copy().replace(np.nan,0)
    
    # Check if columns are in id_mapping_df
    if protein_col not in id_mapping_df:
        module_output += "###ERROR### Group column for summariaztion not found"
        return data_df, id_mapping_df, module_output
    
    if seq_col not in id_mapping_df:
        module_output += "###ERROR### Sequence column for summariaztion not found"
        return data_df, id_mapping_df, module_output
    
    # add new column to id_mapping_df
    id_mapping_df["Summarisation"] = 0
        
    pgs_list = id_mapping_df[protein_col].unique() 
    id_mapping_dict = {pg:list(id_mapping_df[id_mapping_df[protein_col]==pg][seq_col]) for pg in pgs_list}
    
    id_mapping_df.set_index(seq_col,inplace=True)
    
    #for every protein group
    for pg in pgs_list:
        
        precursors_list = [precursor for precursor in id_mapping_dict[pg] if precursor in indices]
        
        if precursors_list:
            temp_df = data_df.loc[precursors_list]
            
            ##on weighted harmonic mean between normalized abundance, correlation and amount of measured samples
            hmean_df=pd.DataFrame(index=temp_df.index,columns=["norm. Abundance","norm. Correlation","norm. Imputation"])
            sums=temp_df.sum(axis=1)
            corrs=temp_df.T.corr()
            corrs_hmean=[weighted_harmonic_mean(corrs.loc[x]) for x in corrs.index]
            temp_mask=filter_mask_z.loc[temp_df.index].map(lambda x: int(str(float(x)).split(".")[-1]))
            measured_precs=(temp_mask==0).sum(axis=1)
            
            hmean_df["norm. Abundance"]=sums/sums.max() #division by zero if peptide not abundant over all fractions => should be excluded in remove_peptides_with_excessing_nans
            hmean_df["norm. Correlation"]=corrs_hmean/np.max(corrs_hmean) #division by zero if all precursors correlate zero
            hmean_df["norm. Imputation"]=measured_precs/measured_precs.max() #division by zero if all values are imputed for precursor => should be excluded in remove_peptides_with_excessing_nans
            hmean_df["hmean"]=[weighted_harmonic_mean(hmean_df.loc[x],weights) for x in hmean_df.index]
        
            #set indices for filter_mask
            if pooling=="mean":
                mask_indices=temp_df.loc[hmean_df.sort_values(by="hmean",ascending=False).index[:top_n_peptides]].index
            else: 
                prec_df=temp_df.loc[hmean_df.sort_values(by="hmean",ascending=False).index[:top_n_peptides]]
                mask_indices=prec_df.index
                           
            if len(temp_df)>=min_peptides:
                id_mapping_df.loc[mask_indices,"Summarisation"] = 1
                if pooling=="mean":
                    pg_profiles_df.loc[pg]=temp_df.loc[hmean_df.sort_values(by="hmean",ascending=False).index[:top_n_peptides]].mean()
                elif pooling=="hmean":
                    pg_profiles_df.loc[pg]=[weighted_harmonic_mean(prec_df[c].values.reshape(len(prec_df))) for c in prec_df.columns]
                elif pooling=="gmean":
                    pg_profiles_df.loc[pg]=[geometric_mean(prec_df[c].values.reshape(len(prec_df))) for c in prec_df.columns]
                #filter_mask.loc[mask_indices,data_df.columns]="s"
            else:
                id_mapping_df.loc[mask_indices,"Summarisation"] = -1

    pg_profiles_df.index.name="Protein.Group"

    return(pg_profiles_df,id_mapping_df,module_output)

def transform_profiles_to_longformat_table(profiles, id_mapping_df, protein_col, seq_col):
    """
    Reshapes wide format profiles with IDs as rows, samples as columns and abundances as values to long-format with ID, sample and values as columns.
    
    author: BZ
    
    Parameters
    ----------
    profiles : pd.DataFrame
        rows: precursors, columns: samples
    id_mapping_df: pd.DataFrame
        information on mapping from precursor to protein (group) level. rows denote individual precursors
    seq_col : str
        name of index level containing precursor Ids. The default is 'Precursor.Id'.
    protein_col : str
        name of index level containing protein Ids. The default is 'Protein.Group'.

    Returns
    -------
    lf_df : pd.DataFrame
        columns: [ID, sample, value]
    id_mapping_df_c: pd.DataFrame
        information on mapping from precursor to protein (group) level. rows denote individual precursors. Index set to seq_col.

    """
    '''
    ap_df=pd.DataFrame(columns=["precursor","protein_group","filename","prec_intensity"])
    ap_df["precursor"]=profiles.index.get_level_values("Precursor.Id").tolist()*len(profiles.columns)
    ap_df["protein_group"]=profiles.index.get_level_values("Protein.Group").tolist()*len(profiles.columns)
    ap_df["filename"]=np.repeat(profiles.columns.get_level_values("sample").tolist(),len(profiles))
    ap_df["prec_intensity"]=profiles.melt()["value"]
    '''
    
    #profiles.reset_index(inplace=True)
    lf_df=profiles.melt(ignore_index=False)
    lf_df.reset_index(inplace=True)
    lf_df.columns = [seq_col,"sample","value"]
    id_mapping_df_c = id_mapping_df.copy().set_index(seq_col)
    ### Add protein_col
    lf_df[protein_col] = [id_mapping_df_c.loc[lf_df.loc[index,seq_col],protein_col] for index in lf_df.index]
    return(lf_df, id_mapping_df_c)

def precursor_sum(data_df, id_mapping_df, module_output, seq_col="Precursor.Id", protein_col="Protein.Group"):
    """
    sums precursor intensities to obtain protein quantities.

    author: BZ/OL
    
    Parameters
    ----------
    data_df : pandas.DataFrame
        rows: precursors, columns: samples
    id_mapping_df: pd.DataFrame
        information on mapping from precursor to protein (group) level. rows denote individual precursors
    module_output: str
        log string to be filled.
    seq_col : str
        name of index level containing precursor Ids. The default is 'Precursor.Id'.
    protein_col : str
        name of index level containing protein Ids. The default is 'Protein.Group'.

    Returns
    -------
    data_prot_df : pd.DataFrame
        rows equal protein groups, columns correspond to samples and values in the dataframe to intensities.
    module_output: str
        log string to be filled.

    """
    id_mapping_df_c = id_mapping_df.copy().set_index(seq_col)
    data_df.index = [id_mapping_df_c.loc[index,protein_col] for index in data_df.index]
    
    #data_df = 2**data_df
    data_df = data_df.groupby(by="index",level=0).sum()
    #data_df = np.log2(data_df).replace([np.inf, -np.inf], np.nan)
    
    module_output += "\n"+time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tProteins summarized to Precursors as sum of Precursor intensities.\n"
    
    return(data_df, module_output)

if __name__ == "__main__":

    module = "SUMMARISATION"
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
    
    data_df_fp, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, mask_description_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = module, 
                                                                                                                                                                           get_metadata=False, 
                                                                                                                                                                           get_id_mapping=True, 
                                                                                                                                                                           get_mask=True, 
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
    
    helper.logger(log_file,time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tLoading input file " + data_df_fp + ".\n")
    try:
        data_df = helper.load_data(data_df_fp)
        helper.logger(log_file,time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tFile loaded successfully.\n")
    except:
        data_df = None
    
    if data_df is None:
        helper.logger(log_file,"###ERROR###: Input data file could not be loaded")
        if exit_module:
            sys.exit("###ERROR###: Input data file could not be loaded")
    
    #if parameters_dict["summarization_method"] in ["mean","gmean","hmean"]:
    #    masked_df = helper.load_mask(data_df_fp)
    #    mask_description_df = helper.load_mask_description(data_df_fp)
    
    # Check if loading was successful, else pass error to log and exit   
    helper.logger(log_file,module_output)
    module_output = ""
    #if metadata_df is None:
    #    sys.exit("###ERROR###: Metadata loading failed.")
    #if "###WARNING###" in module_output:
    #    helper.logger(log_file,"###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])
    #    sys.exit("###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])

    # unlog data
    data_df=helper.unlog_transform(data_df)
    helper.save_data(data_df, out_dir_fp, outname[:-4]+"_precursors_unlog.tsv")
    helper.logger(log_file,time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tUnlogged precursor data saved to " + os.path.join(out_dir_fp, outname[:-4]+"_precursors_unlog.tsv")+".\n")

    try:
        if parameters_dict["summarization_method"] in ["mean","gmean","hmean"]:
            data_processed_df, id_mapping_df, module_output = topN_summarization(data_df, id_mapping_df, masked_df, module_output, 
                                                                                 top_n_peptides=parameters_dict["top_n_peptides"],
                                                                                 weights=parameters_dict["weights"],
                                                                                 pooling=parameters_dict["summarization_method"],
                                                                                 min_peptides=parameters_dict["min_peptides"], 
                                                                                 seq_col=parameters_dict["seq_col"], 
                                                                                 protein_col = parameters_dict["protein_col"])
        elif parameters_dict["summarization_method"]=="sum":
            data_processed_df, module_output = precursor_sum(data_df, id_mapping_df, module_output, 
                                                             seq_col=parameters_dict["seq_col"], 
                                                             protein_col = parameters_dict["protein_col"])
        else:
            data_processed_df = data_df
            module_output += f"###WARNING###: Unknown summarization method: {parameters_dict['summarization_method']} — skipping.\n"
    except Exception as e:
        e_type, e_obj, e_tb = sys.exc_info()
        script_name = os.path.split(e_tb.tb_frame.f_code.co_filename)[1]
        helper.logger(log_file, "###ERROR### Module " +script_name + " ran into an error at line : " + str(e_tb.tb_lineno) + ": " + str(e) + "\n")
        if exit_module:
            sys.exit("Summarisation exited with an error.\n")
    #except ValueError:
    #    raise(ValueError)

    if parameters_dict["summarization_method"] in ["mean","gmean","hmean"]:
        id_mapping_df.to_csv(os.path.join(out_dir_fp,"id_mapping.tsv"),sep="\t")
    
    # Save without log transformation
    helper.save_data(data_processed_df, out_dir_fp, outname[:-4]+"_unlog.tsv")
    module_output += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tUnlogged protein data saved to " + os.path.join(out_dir_fp, outname[:-4]+"_unlog.tsv")+".\n"
    data_processed_df = helper.log_transform(data_processed_df)

    #if parameters_dict["summarization_method"] in ["mean","gmean","hmean"]:
    helper.close_module(data_processed_df, masked_df, mask_description_df, log_file, module_output, out_dir_fp,outname, mask=True)
    #else:
    #    helper.close_module(data_processed_df, masked_df, mask_description_df, log_file, module_output, out_dir_fp,outname, mask=False)
