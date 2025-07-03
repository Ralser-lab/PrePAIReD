"""
Input: wide-format data
Output: Corrected wide-format data


- batchcorrection_combat: applies combat batch correction
- batchcorrection_harmony: applies harmony batch correction
- correct_for_confounders: drift correction per id for a given set of confounders which are statistically significant correlated to the id
- batchwise_drift_correction: applies drift correction per batch for a set of covariates

"""
# %%

#general
import os
import sys
import warnings
import argparse 
import json
import pandas as pd
import numpy as np
import time
from tqdm import tqdm

warnings.filterwarnings("ignore")

#own scripts
# from Python.Snippets import global_variables as global_variables # kept as place holder for logging, will be removed during programming of main 
# from modules_python import module_helper

#visualization
import matplotlib.pyplot as plt

#statistics and modelling
import pycombat #CoAxLab
# import scanpy as sc
from harmony import harmonize
from sklearn.preprocessing import StandardScaler
from scipy import stats
from sklearn import linear_model
from skmisc import loess
from sklearn import metrics


os.chdir(r"c:\\Users\\cleme\\OneDrive - Charité - Universitätsmedizin Berlin\\Projects\\11_PREPROCESSING\\Preprocessing")

##############################################################################
# Batch correction using parametric ComBat
# Author: SB
##############################################################################

def prepare_input_data(raw_data, log_transform=True):
    """
    Prepares data by optionally applying log2 transformation.

    author: SB

    Parameters
    ----------
    raw_data : pd.DataFrame
        Raw data for batch correction.
    log_transform : bool, optional
        Applies log2 transform if True (default True).

    Returns
    -------
    pd.DataFrame
        Prepared data.
    """
    return np.log2(raw_data) if log_transform else raw_data.copy()

def load_batch_metadata(input_data, meta_path, batchID):
    """
    Loads and merges batch metadata aligned with input data samples.

    author: SB

    Parameters
    ----------
    input_data : pd.DataFrame
        DataFrame of input samples.
    meta_path : str
        Path to metadata Excel file.
    batchID : str
        Column identifier for batch information.

    Returns
    -------
    pd.DataFrame
        DataFrame of merged batch metadata.
    """
    metadata = pd.read_excel(meta_path, index_col=0)
    batchdata = metadata[[batchID]].drop_duplicates()
    input_batchdata = pd.DataFrame(index=input_data.index)
    batchdata_merged = pd.merge(input_batchdata, batchdata, how='inner', left_index=True, right_index=True)
    return batchdata_merged

def parametric_combat(input_data, batchdata, batchID):
    """
    Performs batch correction using ComBat method from Scanpy.

    author: SB

    Parameters
    ----------
    input_data : pd.DataFrame
        Input data to batch correct.
    batchdata : pd.DataFrame
        Batch information aligned to samples.
    batchID : str
        Identifier for batch information.

    Returns
    -------
    pd.DataFrame
        Batch-corrected data.
    """
    adata = sc.AnnData(X=input_data.values, 
                       obs=batchdata, 
                       var=pd.DataFrame(index=input_data.columns))

    sc.pp.combat(adata, batchID)

    corrected_data = pd.DataFrame(adata.X, 
                                  index=adata.obs_names, 
                                  columns=input_data.columns,
                                  dtype=float)
    return corrected_data

def combat_batch_correction(data_in, metadata, batchID):
    """
    Integrated pipeline for preparing data, loading metadata, and performing ComBat batch correction.

    author: SB

    Parameters
    ----------
    data_in : pd.DataFrame
        Data requiring batch correction.
    metadata : str
        Path to metadata Excel file.
    batchID : str
        Batch identifier column in metadata.

    Returns
    -------
    pd.DataFrame
        Batch-corrected dataset.
    """

    # Prepare input data
    input_data = prepare_input_data(data_in, log_transform=True)

    # Load metadata and merge
    batchdata = load_batch_metadata(input_data, meta_path=metadata, batchID=batchID)

    # Perform batch correction
    corrected_data = parametric_combat(input_data, batchdata, batchID=batchID)

    return corrected_data


############################################

def correct_for_confounders(profiles_uncorr,metadata_df,datatype_df,sample_col=None,alpha_significance=0.05,model="loess",random_state=None,fitsamples=None):
    """
    drift correction per id for a given set of confounders which are statistically significant correlated to the id. 
    creates a linear model for every protein over every metadata instance, which is significantly associated to the given protein and corrects its abundance by the model
    the W in LOWESS stands for linear, whereas LOESS uses quadratic regression
    NOTE: no multiple testing correction included.
        
    author: BZ
    
    Parameters
    ----------
    profiles_uncorr : pd.DataFrame
        DataFrame with IDs as index, samples as columns and abundances as values.
    metadata_df : pd.DataFrame
        DataFrame with samples as rows and only covariates for correction as columns.
    datatype_df : pd.DataFrame
        DataFrame with covariates (possibly more metadata factors) as rows and at least "DataType" as column with values in ["ordinal","continuous","nominal"]; depending on covariate DataType the applicable statistical test is selected.
    sample_col : str or None, optional
        If profiles_uncorr.columns is a MultiIndex, sample_col should specify the level containing the sample names as in metadata_df.index. The default is None for columns in profiles_uncorr corresponding to rows in metadata_df.
    alpha_significance : float, optional
        significance level for rejection of the null hypothesis in statistical testing. The default is 0.05.
    model : str, optional
        name of the model used for drift correction, one of ["OLS","ridge","lasso","elasticnet","huber","bayesian","loess"]. The default is "loess".
    random_state : int, RandomState instance or None, optional
        Controls ndom stat. The default is None.
    fitsamples : iterable or None, optional
        sample names to use for construction of the model for correction, e.g. only QC samples. The default is None, using all samples.

    Returns
    -------
    profiles: pd.DataFrame
        corrected profiles, same format as profiles_uncorr.
    r2_df: pd.DataFrame
        containing for every ID the association/correlation, pvalue from test on association/correlation and coefficient of determination R2 from the model.

    """
    if model not in ["OLS","ridge","lasso","elasticnet","huber","bayesian","loess"]:
        print("Model " + model + " not implemented yet, no correction applied.")
        return([None,None]) # or sys.exit(0)
    profiles=profiles_uncorr.copy()
    if metadata_df is None:
        metadata=profiles.columns.to_frame().set_index(sample_col)
    else:
        #some numerical adjustments may be made to the metadata
        metadata=metadata_df.copy()
    #extract sample names from profiles which can be matched with metadata
    if sample_col is not None:
        if fitsamples is not None:
            sample_col_index=profiles.columns.names.index(sample_col)
            fitsamples=[x[sample_col_index] for x in fitsamples]
        profiles_columns_orig=profiles.columns.to_frame(index=False).set_index(sample_col,drop=False)
        profiles.columns=profiles.columns.get_level_values(sample_col)
    samples=profiles.columns
    if len(set(metadata.dropna().index) & set(profiles.columns)) < len(profiles.columns):
        print("WARNING: not for all samples in profiles there is metadata without missing values")
    metadata=metadata.loc[samples]
    if fitsamples is not None:
        samples=list(set(fitsamples) & set(samples))
    #collect performances
    r2_df=pd.DataFrame(index=profiles.index,columns=["R2"])
    r2_df.sort_index(inplace=True)
    #iterate over all proteins
    for prot in profiles.index:
        #identify significantly associated metadata variables
        metadata_vars=list()
        for metadata_var in metadata.columns:
            var_samples=metadata.loc[samples,metadata_var].dropna().index
            var_samples_columns=var_samples
            samples_test=[profiles.loc[prot,metadata[metadata[metadata_var]==x].index] for x in metadata[metadata_var].unique()]
            #test for significant associations/correlations
            if datatype_df.loc[metadata_var,"DataType"]=="nominal":
                if len(metadata.loc[var_samples,metadata_var].unique())==2:
                    #transfer metadata variable to boolean if it is not yet
                    if metadata.loc[var_samples,metadata_var].dtype not in [bool,float]:
                        metadata.loc[var_samples,metadata_var]=metadata.loc[var_samples]==metadata.loc[var_samples].iloc[0]
                    #Point biserial correlation + one way ANOVA
                    association,pval_corr=stats.pointbiserialr(metadata.loc[var_samples,metadata_var],profiles.loc[prot,var_samples_columns])
                    pval_test=stats.f_oneway(*samples_test).pvalue
                    pval=pval_test
                    if (pval_test<=alpha_significance and pval_corr>alpha_significance) or (pval_test>alpha_significance and pval_corr<=alpha_significance):
                        print("WARNING: p-value from correlation and ANOVA disagree on neglecting null hypothesis")
                    r2_df.loc[prot,"pointbiserialr " + metadata_var]=association
                elif len(metadata.loc[var_samples].unique())<2:
                    #no association possible with one value per variable
                    association=0
                    pval=1
                else:
                    pval=stats.f_oneway(*samples_test).pvalue
            elif datatype_df.loc[metadata_var,"DataType"]=="ordinal":
                #Spearman correlation + one way ANOVA
                association,pval_corr=stats.spearmanr(metadata.loc[var_samples,metadata_var],profiles.loc[prot,var_samples_columns])
                pval_test=stats.f_oneway(*samples_test).pvalue
                pval=pval_test
                if (pval_test<=alpha_significance and pval_corr>alpha_significance) or (pval_test>alpha_significance and pval_corr<=alpha_significance):
                    print("WARNING: p-value from correlation and ANOVA disagree on neglecting null hypothesis")
                r2_df.loc[prot,"spearmanr " + metadata_var]=association
            elif datatype_df.loc[metadata_var,"DataType"]=="continuous":
                #Pearson correlation
                association,pval=stats.pearsonr(metadata.loc[var_samples,metadata_var],profiles.loc[prot,var_samples_columns])
                r2_df.loc[prot,"pearsonr " + metadata_var]=association
            if pval<alpha_significance:
                metadata_vars.append(metadata_var)
            r2_df.loc[prot,"pvalue "+ metadata_var]=pval
        if len(metadata_vars)==0:
            continue
        #transfer categorical covariates
        categorical_features=datatype_df.loc[metadata[metadata_vars].columns,"DataType"].isin(["nominal"]) #["nominal","ordinal"]
        categorical_features=categorical_features[categorical_features].index
        if len(categorical_features)>0:
            C_cat=pd.get_dummies(metadata[categorical_features].astype(str),drop_first=True)
            C=C_cat.join(metadata[metadata_vars].drop(categorical_features,axis=1))
        else:
            C=metadata[metadata_vars]
        C.dropna(inplace=True)
        #order profiles and C in the same manner
        common_samples=list(set(profiles.columns) & set(C.index))
        if len(common_samples)<len(profiles.columns):
            print("WARNING: not all samples in profiles are found/annotated in metadata. Incomplete samples are being dropped.")
        if len(common_samples) < len(C):
            print("WARNING: not all samples in metadata are found in profiles. Incomplete samples are being dropped.")
        C=C.loc[common_samples]
        profiles=profiles[common_samples]
        if fitsamples is not None:
            C_f=C.loc[fitsamples]
        else:
            C_f=C
        #define samples
        sample_columns=C_f.index
        #create a linear model with significantly associated/correlated metadata variables
        ##OLS: basic, but assumes that variables are independent, tends to overfitting
        ##Ridge Regression: handles multicollinearity by introducing bias in estimator
        ##Lasso regression: similar to Ridge, but with variable selection
        ##PLS: for fewer observations than variables or variables are correlated, decrease correlated variables to smaller number of uncorrelated principal components using PCA
        ##HuberRegressor: less outlier sensitive
        #scaling
        scaler=StandardScaler()
        scaler.fit(C_f)
        C_fz=scaler.transform(C_f)
        C_z=scaler.transform(C)
        #C_fz=C_f
        #C_z=C
        alphas=np.logspace(-3,3,7)
        if model=="OLS":
            estimator=linear_model.LinearRegression()
        elif model=="ridge":
            #estimator=linear_model.Ridge(alpha=alpha)
            estimator=linear_model.RidgeCV(alphas=alphas)           
        elif model=="lasso":
            estimator=linear_model.LassoCV(alphas=alphas,random_state=random_state)
        elif model=="elasticnet":
            estimator=linear_model.ElasticNetCV(alphas=alphas,random_state=random_state)
        elif model=="huber":
            estimator=linear_model.HuberRegressor()
        elif model=="bayesian":
            estimator=linear_model.BayesianRidge()
            #estimator=linear_model.ARDRegression()
        #elif model=="logistic":
        #    estimator=linear_model.LogisticRegression()
        elif model=="GLM":
            pass
        mean_abundance=profiles.loc[prot].mean()
        if model in ["OLS","ridge","lasso","elasticnet","huber","bayesian"]:
            estimator.fit(C_fz,profiles.loc[prot,sample_columns])
            r2_df.loc[prot,"R2"]=estimator.score(C_fz,profiles.loc[prot,sample_columns])
            coefs=estimator.coef_
            n=estimator.intercept_  
            #scatter plot prot abundance ( metadata variable )
            ##plt.figure()
            ##plt.scatter(C_z,profiles.loc[prot,sample_columns])
            #correction
            if fitsamples is not None:
                for d in range(len(coefs)):
                    profiles.loc[prot,fitsamples]=profiles.loc[prot,fitsamples]-coefs[d]*C_z[:,d]
                profiles.loc[prot,fitsamples]=profiles.loc[prot,fitsamples]-n+mean_abundance
            else:
                for d in range(len(coefs)):
                    profiles.loc[prot]=profiles.loc[prot]-coefs[d]*C_z[:,d]
                profiles.loc[prot]=profiles.loc[prot]-n+mean_abundance
        elif model=="loess":
            estimator=loess.loess(C_fz,profiles.loc[prot,sample_columns])
            try:
                estimator.fit()
            except ValueError:
                print("WARNING: Fitting did not work (probably too similar outcomes), protein skipped.")
                continue
            #print(estimator.input_summary())
            #print(estimator.output_summary())
            pred_values=estimator.predict(C_fz).values
            if fitsamples is not None:
                r2_df.loc[prot,"R2"]=metrics.r2_score(profiles.loc[prot,fitsamples],pred_values)
                profiles.loc[prot,fitsamples]=profiles.loc[prot,fitsamples]-pred_values+mean_abundance
            else:
                r2_df.loc[prot,"R2"]=metrics.r2_score(profiles.loc[prot],pred_values)
                profiles.loc[prot]=profiles.loc[prot]-pred_values+mean_abundance
            #scatter plot prot abundance ( metadata variable )
            '''
            plt.figure()
            #plt.scatter(C_fz[:,1],profiles.loc[prot,sample_columns])
            plt.scatter(C_fz[:,0],profiles.loc[prot,sample_columns])
            
            '''
    #reset multitindex columns
    if sample_col is not None:
        profiles.columns=pd.MultiIndex.from_frame(profiles_columns_orig.loc[profiles.columns])
    return([profiles, r2_df])

def batchwise_drift_correction(logprofiles, metadata_df, metadata_metadata, model="loess", covariates=None, sample_col="sample", batch_level="msbatch", fitsamples=None, correct_one_by_one=True, random_state=None):
    """
    apply drift correction per batch
    
    author: BZ
    
    Parameters
    ----------
    logprofiles : pd.DataFrame
        DataFrame with IDs as rows, samples as columns and intensities as values.
    metadata_df : pd.DataFrame
        DataFrame with samples as rows and columns containing covariates for drift correction in columns.
    metadata_metadata : pd.DataFrame
        DataFrame with covariates (possibly more metadata factors) as rows and at least "DataType" as column with values in ["ordinal","continuous","nominal"]; depending on covariate DataType the applicable statistical test is selected.
    model : str, optional
        String specifying model used for drift correction. The default is "loess".
    covariates : iterable, optional
        Iterable, e.g. list specifying covariates which should be corrected by drift. The default is None, corresponding to all columns in metadata_df.
    sample_col : str, optional
        String specifying the name of the metadata_df column containing sample names as in logprofiles.columns. The default is "sample".
    batch_level : str, optional
        column level name with batch information in logprofiles. The default is "msbatch".
    fitsamples : iterable or None, optional
        samples corresponding to a subset of logprofiles.columns with samples used to fit the model for correction. The default is None, using all samples.
    correct_one_by_one : bool, optional
        whether to correct all covariates in one model or coonstruct a model per covariate, which also depends on the model used. The default is True.
    random_state : int, RandomState instance or None, optional
        Controls ndom stat. The default is None.

    Returns
    -------
    logprofiles_corr : pd.DataFrame
        DataFrame with IDs as rows, samples as columns and intensities as values.
    r2_df_final : pd.DataFrame
        DataFrame containing for every ID in logprofiles the model evaluation and statistical test results. 

    """
    #check if covariates are in metadata_metadata
    if len(set(covariates)-set(metadata_metadata.index))>0:
        print(str(set(covariates)-set(metadata_metadata.index)) +" not in metadata_metadata")
        return
    if covariates is None:
        covariates=metadata_df.columns
    logprofiles_corr = logprofiles.copy()
    #assign metadata_df to logprofiles
    if fitsamples is not None: #necessary, because logprofiles_corr.columns could get new levels in next if-statement
        fitsamples_indices=logprofiles_corr.columns.isin(fitsamples)
    if len(set(covariates) - set(logprofiles.columns.names))>0:
        logprofiles_corr=module_helper.annotate_metadata_to_profiles(logprofiles_corr, metadata_df[list(set(metadata_df.columns)-set(logprofiles_corr.columns.names) | set([sample_col]))], sample_col=sample_col)
    #check if covariates in logprofiles_columns
    if len(set(covariates) - set(logprofiles.columns.names))>0:
        print(str(set(covariates) - set(logprofiles.columns.names)) +" neither in metadata nor in logprofiles.columns.names")
        return
    if fitsamples is not None:
        fitsamples=logprofiles_corr.columns[fitsamples_indices]
    r2_dfs=list()
    fitsamples_batch=None
    if correct_one_by_one:
        for covariate in covariates:
            #plot the distribution of the covariate in fitsamples vs rest
            all_values = logprofiles_corr.columns.to_frame(index=False).set_index(sample_col)[covariate]
            fitsample_values = logprofiles_corr[fitsamples].columns.to_frame(index=False).set_index(sample_col)[covariate]
            f,ax=plt.subplots()
            ax.hist(all_values,label="All",density=True,alpha=0.5)
            ax.hist(fitsample_values,label="Fitsamples",density=True,alpha=0.5)
            ax.legend()
            ax.set_xlabel(covariate)
            ax.set_ylabel("Density")
            #f.savefig()
            #ax.boxplot([all_values,fitsample_values],tick_labels=["All","Fitsamples"])
            corrected_logprofiles=list()
            batch_index=list(fitsamples.names).index(batch_level)
            for batch in tqdm(logprofiles_corr.columns.get_level_values(batch_level).unique()):
                uncorr_logprofiles=logprofiles_corr.xs(batch,level=batch_level,axis=1,drop_level=False)
                if fitsamples is not None:
                    #fitsamples=[x[0] for x in filter(lambda x: x[1]=="SP.QC", logprofiles_corr.columns)]
                    #fitsamples_batch=list(set(fitsamples) & set(uncorr_logprofiles.columns))
                    fitsamples_batch = list(filter(lambda x: x[batch_index]==batch,fitsamples))
                    if len(fitsamples_batch)==0:
                        print("\nNo samples for model-fitting present in batch "+str(batch)+"\nPlease adjust fitsamples argument.") 
                        continue
                    elif len(fitsamples_batch)<10:
                        print("\nWARNING: Only "+str(len(fitsamples_batch)) + " samples used for fitting correction model in batch: " + str(batch) + ", covariate: " + str(covariate))
                corr_logprofiles, r2_df = correct_for_confounders(profiles_uncorr=uncorr_logprofiles,metadata_df=pd.DataFrame(uncorr_logprofiles.columns.to_frame(index=False).set_index(sample_col)[covariate]),datatype_df=metadata_metadata,sample_col=sample_col,model=model,random_state=random_state,fitsamples=fitsamples_batch)
                r2_df.rename(columns={"R2":"R2 "+covariate})
                r2_df.columns = [x + " batch"+str(batch) for x in r2_df.columns]
                r2_dfs.append(r2_df)
                corrected_logprofiles.append(corr_logprofiles)
            #concatenate all corrected profiles and apply batch correction
            logprofiles_corr=pd.concat(corrected_logprofiles,axis=1)
    else: #correct for all covariates at once
        corrected_logprofiles=list()
        for batch in tqdm(logprofiles_corr.columns.get_level_values(batch_level).unique()):
            uncorr_logprofiles=logprofiles_corr.xs(batch,level="batch_level",axis=1,drop_level=False)
            if fitsamples is not None:
                #fitsamples=[x[0] for x in filter(lambda x: x[1]=="SP.QC", logprofiles_corr.columns)]
                fitsamples_batch=list(set(fitsamples) & set(uncorr_logprofiles.columns))
                if len(fitsamples_batch)==0:
                    print("\nNo samples for model-fitting present in batch "+str(batch)+"\nPlease adjust fitsamples argument.") 
                    continue 
            corr_logprofiles, r2_df = correct_for_confounders(profiles_uncorr=uncorr_logprofiles,metadata_df=pd.DataFrame(uncorr_logprofiles.columns.to_frame(index=False).set_index("sample")[covariates]),datatype_df=metadata_metadata,sample_col="sample",model=model,random_state=42,fitsamples=fitsamples_batch)
            r2_df.columns = [x + " batch"+str(batch) for x in r2_df.columns]
            r2_dfs.append(r2_df)
            #r2_df.to_csv(qc_folder + savefilename_precursors + appendix + "reprep_normalized_imputed,imputed_loesscorrected_on_all_for_SPBatch_"+str(int(spbatch))+".tsv",sep="\t")
            corrected_logprofiles.append(corr_logprofiles)
            #concatenate all corrected profiles
            logprofiles_corr=pd.concat(corrected_logprofiles,axis=1)
    r2_df_final=pd.concat(r2_dfs,axis=1)
    return([logprofiles_corr, r2_df_final])
    
def batchcorrect_harmony(
            data_df:pd.DataFrame, meta_df:pd.DataFrame, 
            batch_key:list, random_seed:int, **kwargs): 

    """
    Applies the HARMONY batch correction algorithm to the input data.
    The function does not accept Na's.

    author: CD

    As implemented here: https://pypi.org/project/harmony-pytorch/

    This function takes an expression matrix and associated metadata, 
    and performs batch correction using the Harmony algorithm based on the 
    specified batch variables.

    Parameters
    ----------
        data_df (pd.DataFrame) : The data matrix (e.g., gene expression), with rows as samples and columns as features.
        meta_df (pd.DataFrame) : Metadata for the samples, must align with `data_df` by index.
        batch_key (list) : List of column names in `meta_df` that define batch groupings.
        random_seed (int) : Random seed for reproducibility.
        **kwargs: Additional keyword arguments passed to the `harmonize` function

    Returns:
    -------
        Tuple[pd.DataFrame, list]: A tuple containing:
            - `data_processed_df` (pd.DataFrame): Batch-corrected data.
            - `log_list` (list of str): Log messages detailing the correction process.
    """
    log_list = []
    log_list.append(time.strftime("%Y-%m-%d%t%H:%M:%S: "+"Start HARMONY batch correction"))

    data_arr = data_df.to_numpy()
    # Apply harmony algorithm
    data_arr = data_df.to_numpy()
    data_processed_arr = harmonize(
        data_arr, meta_df, batch_key, random_state=random_seed,
        **kwargs)
    data_processed_df = pd.DataFrame(
        data_processed_arr, 
        index=data_df.index, columns=data_df.columns)

    log_list.append(
        time.strftime("%Y-%m-%d%t%H:%M:%S: "\
                      +"Finished HARMONY batch correction"))
    log_string = "\n".join(log_list)

    return data_processed_df, log_string

# %%

## remove before pushing!
data_df_fp = "data/Test_data/Reduced_set/Test_data_reduced.tsv"
param_fp = ""
out_dir_fp = ""
meta_df_fp = "data/Test_data/Additional_Information/Metadata.tsv"

parameters_dict = {
    "batchcorrect_method" : "combat" , # string, options "harmony"
    "batch_cols_meta" : ["MS.Batch"], # list os strings for harmony, combat only allows one
    "random_seed" : 42,
    # "harmony_n_clusters": None,
    # "harmony_max_iter_harmony": 10,
    # "harmony_max_iter_clustering": 200,
    # "harmony_tol_harmony": 1e-4,
    # "harmony_tol_clustering": 1e-5,
    # "harmony_ridge_lambda": 1.0,
    # "harmony_sigma": 0.1,
    # "harmony_n_jobs": -1,
    # "harmony_verbose": False,
    # "harmony_theta": 2.0,
    # "harmony_tau": 0,
    "combat_fitsamples_meta" : ["Sample.Type", "QC-Pool"], # column in meta data and group(s) or None
    "combat_cols_preserve_meta" : ["Batch.Injection"], # column in meta file that should not be touched during batch correction
    "combat_cols_adjust_meta" : ["Injection"], # column in meta file that should not be touched during batch correction
}

# %%

if __name__ == "__main__":

    # Return string for logging from the function
    # for large functions log begfore function state
    
    # # Get data from wrapper by parsing flags
    # parser = argparse.ArgumentParser()
    # parser.add_argument(
    #     "-d", "--datapath", required=True, 
    #     help="Filepath to datafile that will be processed") 
    # parser.add_argument(
    #     "-p", "--parampath", required=True, 
    #     help="Filepath to parameter json") 
    # parser.add_argument(
    #     "-o", "--outpath", required=True, 
    #     help="Filepath for outputfile") 
    # parser.add_argument(
    #     "-md", "--metapath", required=True, 
    #     help="Filepath to metadatafile") 
    # parser.add_argument(
    #     "-ma", "--maskpath", required=True, 
    #     help="Filepath to maskfile") 
    # parser.add_argument(
    #     "-l", "--log_file", required=True, 
    #     help="Filepath to logfile") 
    # parser.add_argument(
    #     "-c", "--counter", required=True, 
    #     help="Module counter") 

    # args = parser.parse_args()

    # data_df_fp = args.datapath
    # param_fp = args.parampath
    # out_dir_fp = args.outpath
    # meta_df_fp = args.metapath
    # mask_df_fp = args.maskpath
    # log_file = args.log_file
    # counter = args.counter
    counter = 0

    outname_str = "data_batchcorrected"

    log_string = time.strftime("%Y-%m-%d%t%H:%M:%S: "+ f"Start batch correction {counter}!\n")
    
    
    # Load data
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: "+"Load precusor data!\n")
    data_df = pd.read_csv(data_df_fp,sep="\t",index_col=0)
    data_df = data_df.reset_index()[["File.Name", "Protein.Group", "PG.Quantity"]].drop_duplicates() # remove in final version
    data_df = data_df.pivot(index="File.Name", columns="Protein.Group", values="PG.Quantity")
    

    # remove before 
    # data_df = data_df.fillna(0)
    print(data_df.shape)
    # 

    data_df = data_df.loc[:, data_df.isna().mean()<0.99]
    print(data_df.shape)
    # data_df = data_df.loc[:, data_df.var()>1]
    data_df = data_df.fillna(0)
    

    # drop wenn pushing
    # from sklearn.impute import KNNImputer
    # data_df.iloc[:,:] = KNNImputer().fit_transform(data_df)
    
    # remove in final version
    # Load metadata
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: "+"Load meta data!\n")
    meta_df = pd.read_csv(meta_df_fp,sep="\t",index_col=0)
    # Load parameter file
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: "+"Load parameter file!\n")
    # with open('Test_folder/4_IMPUTATION/Module_Parameters.json') as f:
    #     parameters_dict = json.load(f)
    # Load Mask
    # log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: "+"Load mask!\n")
    # mask = pd.read_csv(mask_df_fp,sep="\t",index_col=0)

    # -------------------------------------------------------------------------#
    # Check Compatibility of feature and meta data                             #
    # -------------------------------------------------------------------------#
    
    # check if batchcolumns in metafile
    for batch_col in parameters_dict["batch_cols_meta"]:
        if batch_col not in meta_df.columns:
            raise ValueError(f"###ERROR###: Batch Column {batch_col} not in metafile columns!")

    # remove before pushing
    meta_df = meta_df[meta_df["MS.Batch"]!=8]
    meta_df["Injection"] = meta_df["Injection"].sample(frac=1)

    # check if data and meta file consist of same samples
    # only keep intersect on index
    common_idx_list = list(set(data_df.index) & set(meta_df.index))
    n_samples_data_int = data_df.shape[0]
    n_samples_meta_int = meta_df.shape[0]
    drop_samples_data_lst = list(set(data_df.index) - set(meta_df.index))
    drop_samples_meta_lst = list(set(meta_df.index) - set(data_df.index))

    data_df = data_df.loc[common_idx_list]
    meta_df = meta_df.loc[common_idx_list]
    
    flag_sample_intersect_int = False
    warning_intersect_data_str = ""
    if n_samples_data_int != data_df.shape[0]:
        flag_sample_intersect_int = True
        n_final_samples_data_int = n_samples_data_int - data_df.shape[0]
        warning_intersect_samples_str = \
                time.strftime('%Y-%m-%d\t%H:%M:%S: ') \
                +f"###WARNING###: {n_final_samples_data_int} sample(s) did not appear in meta data:"\
                + "\nFollowing sample(s) has/have been dropped for analysis:" \
                + "\n\t"+'\n\t'.join(drop_samples_data_lst) \
                + "\n"
        log_string += warning_intersect_samples_str
    if n_samples_meta_int != meta_df.shape[0]:
        flag_sample_intersect_int = True
        n_final_samples_meta_int = n_samples_meta_int - meta_df.shape[0]
        warning_intersect_samples_str = \
                time.strftime('%Y-%m-%d\t%H:%M:%S: ') \
                +f"###WARNING###: {n_final_samples_meta_int} sample(s) did not appear in precursor data:\n"\
                + "\n\t\t Following sample(s) has/have been dropped for analysis:\n" \
                + "\t\t" + "\n\t\t".join(drop_samples_meta_lst) \
                + "\n"
        log_string += warning_intersect_samples_str
    
    if flag_sample_intersect_int:
        log_string += time.strftime('%Y-%m-%d\t%H:%M:%S ') \
                +f"{meta_df.shape[0]} remain for analysis.\n"
        


    # Warn user that some samples do not have a batch annotation
    for batch_col in parameters_dict["batch_cols_meta"]:
        samples_na_list = []
        if meta_df[batch_col].isna().sum() > 0:
            count_na_int = meta_df[batch_col].isna().sum()
            samples_na_list = list(
                meta_df[meta_df[batch_col].isna()].index)
            warning_na_meta_str = \
                time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                +f"###WARNING###: {count_na_int} sample(s) do(es) not have an {batch_col} value:"\
                + "\nFollowing sample(s) has/have been dropped for analysis:" \
                + "\n\t"+'\n\t'.join(samples_na_list) + "\n"
            log_string += warning_na_meta_str

            # drop samples
            if len(samples_na_list) > 0:
                drop_samples = set(samples_na_list) & set(data_df.index)
                data_df = data_df.drop(drop_samples, axis=0)
                meta_df = meta_df.drop(drop_samples, axis=0)

    # Exit module if there are multiple levels in batch_cols of meta file
    for b in parameters_dict["batch_cols_meta"]:
        if len(np.unique(meta_df[b])) == 1:
            warning_batch_level_str = \
                time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                +f"###ERROR###: Batch Level {b} does only have 1 level:{np.unique(meta_df[b])[0]} value:"\
                + "\t\tMultiple levels are required to perform batch correction"  + "\n"
            log_string += warning_batch_level_str
            print(log_string)
            sys.exit(0)
    
    
    # -------------------------------------------------------------------------#
    # Perform Batch Correction                                                 #
    # -------------------------------------------------------------------------#

    if parameters_dict["batchcorrect_method"] == "harmony":

        # translate harmony parameter from dict 
        harmony_kwargs_dict = {
                "_".join(k.split("_")[1:]):v for k,v in parameters_dict.items()
                if k.startswith("harmony")}

        data_processed_df, log_string_harmony = batchcorrect_harmony(
            data_df, meta_df, 
            parameters_dict["batch_cols_meta"],
            random_seed=parameters_dict["random_seed"] + counter, 
            **harmony_kwargs_dict
        )
        log_string += log_string_harmony

    
    elif parameters_dict["batchcorrect_method"] == "combat":

        # Check if only one batch col is defined!
        if len(parameters_dict["batch_cols_meta"]) > 1:
            info_fitsamples_str = \
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"COMBAT is only implemented for one batch variable.\n"\
                    +f"{len(parameters_dict['batch_cols_meta'])} levels are provided: {parameters_dict['batch_cols_meta']}"  + "\n"
            log_string += info_fitsamples_str
            print(log_string)
            sys.exit(0)

        # translate combat parameter from dict 
        combat_kwargs_dict = {
                "_".join(k.split("_")[1:]):v for k,v in parameters_dict.items()
                if k.startswith("combat")}
        
        fitsamples_idx_list = None
        if len(combat_kwargs_dict["fitsamples_meta"]) == 2:
            fit_col, fit_level = combat_kwargs_dict["fitsamples_meta"]
            fitsamples_idx_list = list(meta_df[meta_df[fit_col]==fit_level].index)
            info_fitsamples_str = \
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"COMBAT will be fitted on {len(fitsamples_idx_list)} samples.\n"\
                    +f"Level {fit_level} in meta column {fit_col} will be used to fit COMBAT Model."  + "\n"
            log_string += info_fitsamples_str

            # Check whether the same batch levels are defined in samples
            # fitting combat and the ones and the rest of the dataset

            meta_df[fit_col] 

            levels_fitsamples = set(
                meta_df.loc[
                    fitsamples_idx_list, parameters_dict["batch_cols_meta"][0]])
            levels_data = set(
                meta_df[parameters_dict["batch_cols_meta"][0]])
            if len(levels_data) != len(levels_fitsamples):
                info_fitsamples_str = \
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###WARNING###: Uneven numbers of batch levels in whole dataset (n={len(levels_data)}) and samples selected to fit COMBAT model (n={len(levels_fitsamples)}).\n"\
                    +f"Level(s) {list(levels_data-levels_fitsamples)} do not appear in samples to fit model and samples will be discarded"  + "\n"
                log_string += info_fitsamples_str
            # Remove samples of batch levels do not appear in fit samples
            meta_df = meta_df[
                ~meta_df[parameters_dict["batch_cols_meta"][0]].isin(
                    list(levels_data-levels_fitsamples))]
            data_df = data_df.loc[meta_df.index]


        def batchcorrect_combat(
            data_df:pd.DataFrame, 
            meta_df:pd.DataFrame, 
            batch_key:str,
            fitsamples_idx_list=list,
            cov_preserve_list=list,
            cov_adjust_list=list,
            ):
                
                """
                Apply combat batch correction, see https://github.com/CoAxLab/pycombat
                
                #TODO
                check whether covariates provided overlap with profiles.columns, 
                if not inform the user
                how to handle fit samples

                To low variation might lead to a dataframe full of NA's
                if so recheck imputation and variance thresholds


                author: BZ
                
                Parameters
                ----------
                profiles : pd.DataFrame
                    features (e.g. precursors) as rows, samples as columns
                batchkey : str
                    column level to correct by (e.g. "batch")
                fitsamples : list, optional
                    samples to use for fitting the combat class (e.g. QC-samples only), 
                    must be column names of profiles. 
                    The default is None, fitting with all samples in profiles.
                covariates : list
                    list of column levels with covariates to be corrected for
                binlevel : bool
                    if True, bin the levels provided by batchkey
                bins : int
                    number of bins for the binning of batchkey
                sample_column : str
                    name of the column level with the sample names. The default is "sample".
                    
                Returns
                -------
                pd.DataFrame, same shape as input profiles.

                """

                ### Function 

                log_list = []
                log_list.append(
                    time.strftime(
                        "%Y-%m-%d%t%H:%M:%S: "+"Start COMBAT batch correction"))

                # effects to preserve during batch correction
                if cov_preserve_list != None:
                    x_df = meta_df[cov_preserve_list]

                # covariates to adjust for during batch correction
                if cov_adjust_list != None:
                    c_df = meta_df[cov_adjust_list]

                # Init combat Model
                combat_model=pycombat.Combat()
                # fit combat model on selected samples
                if fitsamples_idx_list != None:
                    
                    log_list.append(
                        time.strftime(
                            "%Y-%m-%d%t%H:%M:%S: "+"Fit COMBAT Class on defined samples."))

                    c_fit = combat_model.fit(
                        Y=data_df.loc[fitsamples_idx_list].to_numpy(),
                        b=meta_df[batch_key].loc[fitsamples_idx_list].to_numpy(),
                        X=x_df.loc[fitsamples_idx_list], 
                        C=c_df.loc[fitsamples_idx_list]
                    )
                else:
                    
                    log_list.append(
                        time.strftime(
                            "%Y-%m-%d%t%H:%M:%S: "+"Fit COMBAT Class on all samples."))

                    c_fit = combat_model.fit(
                        Y=data_df.to_numpy(),
                        b=meta_df[batch_key].to_numpy(),
                        X=x_df, 
                        C=c_df
                    )

                log_list.append(
                        time.strftime(
                            "%Y-%m-%d%t%H:%M:%S: "+"Transform samples."))

                data_processed_arr = c_fit.transform(
                    Y=data_df.to_numpy(),
                    b=meta_df[batch_key].to_numpy(),
                    X=x_df, 
                    C=c_df
                )
                data_processed_df = pd.DataFrame(
                    data_processed_arr,
                    index=data_df.index,
                    columns=data_df.columns)

                # Check if output df consists of NA's exclusively 
                # if so return error in Log and don't return processed df
                if data_processed_df.isna().all().all():
                    log_list.append(
                        time.strftime(
                            "%Y-%m-%d%t%H:%M:%S: "+"###ERROR###: All Values in DF were transformed to NA.\n"\
                            +f"Please filter input dataframe for features with high missing values or low variance."  + "\n"))
                    log_string = "\n".join(log_list)
                    return log_string
                else:
                    log_list.append(time.strftime(
                        "%Y-%m-%d%t%H:%M:%S: "\
                        +"Finished HARMONY batch correction successfully."))
                    log_string = "\n".join(log_list)
                    return data_processed_df, log_string
        


        out = batchcorrect_combat(
            data_df = data_df,
            meta_df = meta_df,
            batch_key = parameters_dict["batch_cols_meta"][0],
            fitsamples_idx_list=fitsamples_idx_list,
            cov_preserve_list = parameters_dict["combat_cols_preserve_meta"],
            cov_adjust_list = parameters_dict["combat_cols_adjust_meta"],
        )
        if isinstance(out, str) and "ERROR" in out:
            print(log_string)
            sys.exit(0)
            

        
    # elif parameters_dict["normalization_method"] == "median_scaling_intensity":
    #     data_processed_df = normalize_median_scaling_intensity(data_df, parameters_dict["axis"])
        
    # elif parameters_dict["normalization_method"] == "median":
    #     data_processed_df = normalize_median(data_df)
        
    # elif parameters_dict["normalization_method"] == "mean":
    #     data_processed_df = normalize_mean(data_df)
        
    # elif parameters_dict["normalization_method"] == "delayed":
    #     data_processed_arr = normalize_delayed(data_arr, parameters_dict["method"], parameters_dict["bounds"])

    # # Comparison for mask and logger what changed between both files
    # """
    # # Helper DataFrames
    # mask = helper.update_mask(mask)
    # log = helper.logger(logf)
    # """
    
    # # Save data
    # if parameters_dict["normalization_method"] not in ["delayed"]:
    #     data_processed_df.to_csv(os.path.join(out_dir_fp,outname+".tsv"),sep="\t")
    # else:
    #     np.save(os.path.join(out_dir_fp,out_dir_fp+".npy"),data_processed_arr)

    print(log_string)
    # sys.exit(0)
    
    # Final Method to evaluate processed data df

# %%
