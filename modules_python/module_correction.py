"""
Input: wide-format data
Output: Corrected wide-format data


- batchcorrection_combat: applies combat batch correction
- batchcorrection_harmony: applies harmony batch correction
- get_min_samples: calculates the minimum number of samples needed to calculate covariate association and apply correction by a rule of thumb (10*n_covariates)
- get_ids_for_correction: calculates association and p-value of IDs with covariates to obtain relevant IDs for covariate correction
- adjust_confounders: drift correction per id for a given set of confounders which are statistically significant correlated to the id
- adjust_confounders_batchwise: applies drift correction per batch for a set of covariates

General comments:
- No duplicated adjustment in Combat and via linear modeling
- Other rules of thumb for setting minimum number of samples needed to apply correction can be implemented

# TBD save of corr_dfs in adjust_confounders and adjust_confounders_batchwise missing
# TBD Detailed Module Overview
# Implemented for features as columns and samples as rows

"""
# %%

#general
import sys
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import time
import os
from tqdm import tqdm

#statistics and modelling
import pycombat #CoAxLab
from harmony import harmonize
from sklearn.preprocessing import StandardScaler
from scipy import stats
from sklearn import linear_model
from skmisc import loess
from sklearn import metrics

import module_helper as helper

def batchcorrect_harmony(
            data_df:pd.DataFrame, metadata_df:pd.DataFrame, 
            batch_key:list, seed, **kwargs): 

    """
    This function takes an expression matrix and associated metadata, 
    and performs batch correction using the Harmony algorithm based on the 
    specified batch variables.
    The function does not accept Na's.

    author: BZ, CD

    As implemented here: https://pypi.org/project/harmony-pytorch/

    Parameters
    ----------
        data_df (pd.DataFrame) : The data matrix (e.g., gene expression), with rows as samples and columns as features.
        metadata_df (pd.DataFrame) : Metadata for the samples, must align with `data_df` by index.
        batch_key (list) : List of column names in `metadata_df` that define batch groupings.
        seed (int) : Random seed for reproducibility.
        **kwargs: Additional keyword arguments passed to the `harmonize` function

    Returns:
    -------
        Tuple[pd.DataFrame, list]: A tuple containing:
            - `data_processed_df` (pd.DataFrame): Batch-corrected data.
            - `log_list` (list of str): Log messages detailing the correction process.
    """
    log_string = ""
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: "+"\tStart HARMONY batch correction\n")

    data_arr = data_df.to_numpy()
    # Apply harmony algorithm
    data_arr = data_df.to_numpy()
    data_processed_arr = harmonize(
        data_arr, metadata_df, batch_key, random_state=seed,
        **kwargs)
    data_processed_df = pd.DataFrame(
        data_processed_arr, 
        index=data_df.index, columns=data_df.columns)

    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S"\
                      +"\tFinished HARMONY batch correction")

    return data_processed_df, log_string

def batchcorrect_combat(
    data_df:pd.DataFrame, 
    metadata_df:pd.DataFrame, 
    batch_key:str,
    fitsamples_idx_list=list,
    cov_preserve_list=list,
    cov_adjust_list=list,
    cov_adjust_list_dtypes=list,
    ):
        
    """
    Apply ComBat batch correction by using the pyComBat package.
    See https://github.com/CoAxLab/pycombat for more information.
    
    This function supports optional selection of samples for model 
    fittingand allows specification of covariates to preserve or 
    adjust for during correction.

    author: BZ, CD

    Notes:
    -----
    - If the corrected data contains only NaNs, the function returns an error message in the log instead.
    - Ensure that the covariate columns in `cov_preserve_list` and `cov_adjust_list` exist in `metadata_df`.
    - Low variance in features or improper imputation can lead to failed correction with all NaN results.
    
    Parameters
    ----------
    data_df : pd.DataFrame
        features (e.g. precursors) as rows, samples as columns
    batch_key : str
        column level to correct by (e.g. "batch")
    fitsamples_idx_list : list, optional
        List of sample indices to use for fitting the combat class (e.g. QC-samples only), 
        Must be column names of profiles. 
        The default is None, fitting with all samples.
    cov_preserve_list : list
        list of column levels with covariates to be left untouched
    cov_adjust_list : list
        list of column levels with covariates to be corrected for
    cov_adjust_list_dtypes : list 
        list of dtypes for covariates

    Returns
    -------
    Corrected feature matrix as pd.DataFrame, same shape as input profiles.
    """
    log_string = ""
    log_string += time.strftime(
            "%Y-%m-%d%t%H:%M:%S"+"\tStart ComBat batch correction.\n")
    
    missing_values = data_df.isna()
    if missing_values.sum().sum()>0:
        log_string += time.strftime(
                "###WARNING###: Missing values detected for ComBat. Inserting feature median per batch temporary. The imputed values will be removed after this step.\n")
    
    for item in metadata_df[batch_key].unique():
        indices_batch = metadata_df[metadata_df[batch_key]==item].index
        imp_tmp = data_df.loc[indices_batch]
        imp_tmp = imp_tmp.fillna(imp_tmp.median())
        data_df.loc[imp_tmp.index] = imp_tmp
    # If a precursor is missing in full batch set it to median temporarily
    data_df = data_df.fillna(data_df.median())
                
    # effects to preserve during batch correction
    if cov_preserve_list is not None:
        x_df = metadata_df[cov_preserve_list]
    else:
        x_df = None

    # covariates to adjust for during batch correction
    if cov_adjust_list is not None:
        c_df = metadata_df[cov_adjust_list]
        c_vars_nom = [
            c for c, dtype in zip(
                cov_adjust_list, cov_adjust_list_dtypes) if dtype == "nominal"]
        c_vars_not_nom = [
            c for c, dtype in zip(
                cov_adjust_list, cov_adjust_list_dtypes) if dtype != "nominal"]
        
        if  len(c_vars_nom) > 0:
            c_df = pd.concat([
                c_df[c_vars_not_nom], 
                pd.get_dummies(c_df[c_vars_nom].astype(str), drop_first=True)
                        ], axis=1)
    else:
        c_df = None

    # Init combat Model
    combat_model=pycombat.Combat()
    # fit combat model on selected samples
    if fitsamples_idx_list is not None:
        log_string += time.strftime(
                "%Y-%m-%d%t%H:%M:%S"+"\tFit ComBat Class on defined samples.\n")

        Y_fit = data_df.loc[fitsamples_idx_list].to_numpy()
        b_fit = metadata_df[batch_key].loc[fitsamples_idx_list].astype(str).to_numpy()
        if x_df is not None:
            X_fit = x_df.loc[fitsamples_idx_list]
        else:
            X_fit = None
        if c_df is not None:
            C_fit = c_df.loc[fitsamples_idx_list]
        else:
            C_fit = None

        c_fit = combat_model.fit(
            Y=Y_fit,
            b=b_fit,
            X=X_fit,
            C=C_fit
        )
    else:
        log_string += time.strftime(
                "%Y-%m-%d%t%H:%M:%S"+"\tFit ComBat Class on all samples.\n")

        c_fit = combat_model.fit(
            Y=data_df.to_numpy(),
            b=metadata_df[batch_key].to_numpy(),
            X=x_df, 
            C=c_df
        )

    log_string += time.strftime(
                "%Y-%m-%d%t%H:%M:%S"+"\tTransform data.\n")

    samples_all = data_df.index
    Y_all = data_df.to_numpy()
    
    if fitsamples_idx_list is not None:
        b_all =  metadata_df[batch_key].loc[samples_all].astype(str).to_numpy()
    else:
        b_all = metadata_df[batch_key].to_numpy()
    

    X_all = x_df.loc[samples_all] if x_df is not None else None
    C_all = c_df.loc[samples_all] if c_df is not None else None

    # Check if all samples have the same number of observations
    if X_all is not None and C_all is not None:
        assert Y_all.shape[0] == X_all.shape[0] == C_all.shape[0] == b_all.shape[0]
    elif X_all is not None:
        assert Y_all.shape[0] == X_all.shape[0] == b_all.shape[0]
    elif C_all is not None:
        assert Y_all.shape[0] == C_all.shape[0] == b_all.shape[0]
    else:
        assert Y_all.shape[0] == b_all.shape[0]

    data_processed_arr = c_fit.transform(
        Y=Y_all,
        b=b_all,
        X=X_all,
        C=C_all
    )

    data_processed_df = pd.DataFrame(
        data_processed_arr,
        index=data_df.index,
        columns=data_df.columns)
    
    data_processed_df[missing_values] = np.nan

    # Check if output df consists of NA's exclusively
    # if so return error in Log and don't return processed df
    if data_processed_df.isna().all().all():
        log_string += "###ERROR###: All Values in DF were transformed to NA.\n"\
                +"\t\t\t\tPlease filter input dataframe for features with high missing values or low variance."  + "\n"
        return log_string
    else:
        log_string += time.strftime(
            "%Y-%m-%d%t%H:%M:%S"\
            +"\tFinished COMBAT batch correction successfully.\n")

        return data_processed_df, log_string
        
def get_min_samples(m):
    """
    define minimum number of samples needed for fitting 
    rule of thumb: 10 times number of variables, 
    E.g. if fitting a model with 2 variables (InjectionIndex and NumberOfPrecursors simultaneously) 
    at least 20 samples are needed for fitting the model 
    and calculating the association of precursors with these two variables.
    
    #see #https://stats.stackexchange.com/questions/116132/regression-with-very-small-sample-size

    Parameters
    ----------
    m : int
        number of variables in the model.

    Returns
    -------
    int
        minimal number of samples needed to calculate association, statistical test or fit a model.

    """
    return(m * 10)

def get_ids_for_correction(data_df,metadata_df,adjust_cols_list,adjust_cols_dtypes_list,fitsamples_idx_list=None,alpha_significance=0.05):
    
    """
    Calculates association and p-value of IDs (e.g. precursors) with covariates to obtain relevant IDs for confounder correction.
    

    author: BZ, CD

    Parameters:
    -----------
    data_df : pandas.DataFrame
        DataFrame containing the feature data (e.g., gene expression).
        Rows are samples, columns are features.

    metadata_df : pandas.DataFrame
        DataFrame containing metadata variables (e.g., age, sex, batch).
        Must include all columns listed in `adjust_cols_list`.

    adjust_cols_list : list of str
        List of metadata column names to consider for adjustment.

    adjust_cols_dtypes_list : list of str
        List of corresponding data types for each column in `adjust_cols_list`.
        Options: 'nominal', 'ordinal', or 'continuous'.
        
    fitsamples_idx_list : list-like, optional
        Optional list of sample indices to use for model fitting.
        If None, all samples in `data_df` are used.

    alpha_significance : float, optional
        Significance threshold (e.g., 0.05) for determining if a metadata variable is associated 
        with a feature and should be used for adjustment.
        
    Returns:
    --------

    corr_df : pandas.DataFrame
        DataFrame summarizing the association statistics between metadata and features.
        Includes correlation coefficients and p-values for each tested variable.
        If no significant association was found, NaN is reported for r and p-value.
    
    log_string : str
        String summarizing the adjustment process, warnings, and errors encountered.

    Notes:
    ------

    - Incomplete or mismatched sample indices between `data_df` and `metadata_df` are logged and dropped.
    - Uses ANOVA and correlation metrics (point-biserial, Spearman, Pearson) for assessing association:
        - Nominal variables: ANOVA and point-biserial
        - Ordinal variables: ANOVA and Spearman
        - Continuous variables: Pearson

    Warnings:
    ---------
    - Skips features if metadata groups contain fewer than 5 samples.
    - Disagreement between ANOVA and correlation test significance is reported.

    """
    
    log_string = ""
    
    # select samples used for fit model
    if fitsamples_idx_list is None:
        samples_idx_list = data_df.index
    else:
        samples_idx_list = fitsamples_idx_list
    
    # collect associations with meta data variables
    corr_df=pd.DataFrame(index=data_df.columns)
    corr_df.sort_index(inplace=True)

    #define minimum number of samples needed for fitting (rule of thumb)
    min_samples = get_min_samples(len(adjust_cols_list))
            
    # iterate over features
    for idx in data_df.columns:
        metadata_vars = [] # collect variables that show significant associations with feature
        
        
        for metadata_var, dt in zip(adjust_cols_list,adjust_cols_dtypes_list):
            # remove missing values and define common samples in meta and feature data
            samples_meta = metadata_df.loc[samples_idx_list, metadata_var].dropna()
            samples_data = data_df.loc[samples_idx_list, idx].dropna()
            common_samples = list(set(samples_meta.index) & set(samples_data.index))

            samples_idx_test_list= [list(metadata_df[metadata_df[metadata_var]==x].index) for x in metadata_df[metadata_var].unique()]
            samples_idx_test_list = [list(set(x) & set(common_samples)) for x in samples_idx_test_list]
            
            # correlation coefficient to be calculated depending on datatype
            coef_type_str = ""
            if dt == "nominal":
                if len(metadata_df.loc[common_samples,metadata_var].unique())==2:
                    coef_type_str = "pointbiserialr_"
                else:
                    coef_type_str = "anova_"
            elif dt == "ordinal":
                coef_type_str = "spearmanr_"
            elif dt == "continuous":
                coef_type_str = "pearsonr_"
            else:
                pass

            # check how many samples are available per feature
            # do not check association if < min_samples
            if len(common_samples) < min_samples:
                # log_string += f"###WARNING###: Less than 5 samples in data with annotated {metadata_var} in metadata.\n" \
                #         f"\t\t\t\tNo Correction will be applied on {idx}.\n"
                corr_df.loc[idx, "pval_"+ metadata_var] = np.nan
                corr_df.loc[idx, coef_type_str+ metadata_var] = np.nan
                continue


            # check if ordinal or nominal variable groups have at 
            # do not check association if < min_samples samples per group
            if (dt == "ordinal") or (dt=="nominal"):
                if any([len(x) for x in samples_idx_test_list]) < min_samples:
                    log_string += f"###WARNING###: One of the groups in col {metadata_var} consists of less than {min_samples} samples.\n" \
                        f"\t\t\t\tNo Correction will be applied on {idx}.\n"

                    
                    # save nan if sample size was not sufficient
                    corr_df.loc[idx, "pval_"+ metadata_var] = np.nan
                    corr_df.loc[idx, coef_type_str+ metadata_var] = np.nan
                    continue
            
            samples_test_list = [data_df.loc[x, idx] for x in samples_idx_test_list]
            #test for significant associations/correlations
            if dt=="nominal":
                if len(metadata_df.loc[common_samples,metadata_var].unique())==2:
                    # cast metadata variable to boolean for pointbiserial
                    if metadata_df.loc[common_samples,metadata_var].dtype not in [bool, float]:
                        levels_var = metadata_df.loc[common_samples,metadata_var].unique()
                        metadata_df[metadata_var] = metadata_df[metadata_var].replace(
                            {levels_var[0]:True, levels_var[1]:False})
                    # Point biserial correlation + one way ANOVA
                    
                    association, pval_corr = stats.pointbiserialr(
                        metadata_df.loc[common_samples,metadata_var],
                        data_df.loc[common_samples, idx])
                    pval_test = stats.f_oneway(*samples_test_list).pvalue
                    pval = pval_test
                    if (pval_test<=alpha_significance and pval_corr>alpha_significance) or (pval_test>alpha_significance and pval_corr<=alpha_significance):
                        log_string += f"###WARNING###: p-value from correlation and ANOVA disagree on neglecting null hypothesis for {idx}\n"
                    corr_df.loc[idx,coef_type_str + metadata_var] = association
                elif len(metadata_df.loc[common_samples,metadata_var].unique())<2:
                    #no association possible with one value per variable
                    association=0
                    pval=1
                else:
                    pval=stats.f_oneway(*samples_test_list).pvalue
            
            elif dt == "ordinal":
                #Spearman correlation + one way ANOVA
                association, pval_corr = stats.spearmanr(
                    metadata_df.loc[common_samples, metadata_var],
                    data_df.loc[common_samples, idx])
                pval_test=stats.f_oneway(*samples_test_list).pvalue
                pval = pval_test
                if (pval_test<=alpha_significance and pval_corr>alpha_significance) or (pval_test>alpha_significance and pval_corr<=alpha_significance):
                    log_string += f"###WARNING###: p-value from correlation and ANOVA disagree on neglecting null hypothesis for {idx}\n"
                corr_df.loc[idx,coef_type_str + metadata_var]=association

            elif dt =="continuous":
                #Pearson correlation
                association, pval=stats.pearsonr(
                    metadata_df.loc[common_samples, metadata_var],
                    data_df.loc[common_samples, idx])
                corr_df.loc[idx,coef_type_str + metadata_var] = association
            else:
                pass

            # if p-value of association test below significance value 
            # variable will be stored for later adjustment
            if pval < alpha_significance:
                metadata_vars.append(metadata_var)
            corr_df.loc[idx, "pval_"+ metadata_var] = pval

    return corr_df, log_string

def adjust_confounders(
                data_df, metadata_df, model_string, adjust_cols_list,
                adjust_cols_dtypes_list, random_state, fitsamples_idx_list=None,
                alpha_significance=0.05, alphas_models=[-3,3,7], corr_df=None
            ):
        
    """
    Adjusts features in a data matrix by removing the influence of confounding variables 
    from the metadata using a specified regression model. It performs feature-wise analysis 
    and correction based on statistically significant associations between metadata variables 
    and data features.

    author: BZ, CD

    Parameters:
    -----------
    data_df : pandas.DataFrame
        DataFrame containing the feature data (e.g., gene expression).
        Rows are samples, columns are features.

    metadata_df : pandas.DataFrame
        DataFrame containing metadata variables (e.g., age, sex, batch).
        Must include all columns listed in `adjust_cols_list`.

    model_string : str
        Specifies which regression model to use for adjustment. Options:
        - 'ols': Ordinary Least Squares (basic linear regression)
            - Assumes independent variables; can overfit
        - 'ridge': Ridge Regression
            - Adds L2 penalty to reduce multicollinearity
        - 'lasso': Lasso Regression
            - Adds L1 penalty; performs variable selection
        - 'elasticnet': ElasticNet Regression
            - Combines L1 and L2 penalties; balances ridge and lasso
        - 'huber': Huber Regressor
            - Robust to outliers
        - 'bayesian': Bayesian Ridge Regression
            - Incorporates prior distributions; regularizes automatically
        - 'logistic': Logistic Regression
            - Classification-based model (not recommended for continuous feature adjustment)
        - 'lowess': Local polynomial regression (non-parametric, robust to non-linearity)

    adjust_cols_list : list of str
        List of metadata column names to consider for adjustment.

    adjust_cols_dtypes_list : list of str
        List of corresponding data types for each column in `adjust_cols_list`.
        Options: 'nominal', 'ordinal', or 'continuous'.

    random_state : int
        Random state for reproducibility in models that require it (e.g., LassoCV, ElasticNetCV).

    fitsamples_idx_list : list-like, optional
        Optional list of sample indices to use for model fitting.
        If None, all samples in `data_df` are used.

    alpha_significance : float, optional
        Significance threshold (e.g., 0.05) for determining if a metadata variable is associated 
        with a feature and should be used for adjustment.

    alphas_models : list-like, optional
        List of alpha values for models with regularization (e.g., RidgeCV, LassoCV, ElasticNetCV).
    
    corr_df : pandas.DataFrame, optional
        DataFrame summarizing the association statistics between metadata and features.
        Includes correlation coefficients and p-values for each tested variable.
        If no significant association was found, NaN is reported for r and p-value.
        
    Returns:
    --------
    data_processed_df : pandas.DataFrame
        The adjusted feature matrix with confounding effects removed.

    log_string : str
        String summarizing the adjustment process, warnings, and errors encountered.

    corr_df : pandas.DataFrame
        DataFrame summarizing the association statistics between metadata and features.
        Includes correlation coefficients and p-values for each tested variable.
        If no significant association was found, NaN is reported for r and p-value.

    Notes:
    ------
    - Features significantly associated with metadata variables (based on `alpha_significance`) 
      are adjusted using the selected regression method.
    - Categorical metadata variables are one-hot encoded for regression modeling.
    - Model fitting is performed using standardized (z-scored) metadata.
    - Incomplete or mismatched sample indices between `data_df` and `metadata_df` are logged and dropped.
    - Uses ANOVA and correlation metrics (point-biserial, Spearman, Pearson) for assessing association:
        - Nominal variables: ANOVA and point-biserial
        - Ordinal variables: ANOVA and Spearman
        - Continuous variables: Pearson
    
    Warnings:
    ---------
    - Skips features if metadata groups contain fewer than 5 samples.
    - Disagreement between ANOVA and correlation test significance is reported.

    """

    data_processed_df = data_df.copy()

    log_string = ""
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tStart Adjustment.\n"
    
    #calculate the association and p-value for all IDs if not provided as function argument
    if corr_df is None:
        corr_df, log_str = get_ids_for_correction(data_df=data_processed_df,
                                                  metadata_df=metadata_df,
                                                  adjust_cols_list=adjust_cols_list,
                                                  adjust_cols_dtypes_list=adjust_cols_dtypes_list,
                                                  fitsamples_idx_list=fitsamples_idx_list,
                                                  alpha_significance=alpha_significance)
        log_string += log_str
        
    
    # iterate over features
    for idx in data_df.columns:
        metadata_vars = [] # collect variables that show significant associations with feature
        for metadata_var, dt in zip(adjust_cols_list,adjust_cols_dtypes_list):
            if corr_df.loc[idx, "pval_"+ metadata_var] < alpha_significance:
                metadata_vars.append(metadata_var)
        if len(metadata_vars)==0:
            continue
        
        min_samples = get_min_samples(m=len(metadata_vars))
        
        
        #transfer categorical covariates for model fitting
        idx_catgorical_list = [
            i for i, d in enumerate(
                adjust_cols_dtypes_list) if d == "nominal"]
        col_cat_list = [adjust_cols_list[i] for i in idx_catgorical_list]
        col_cat_list =  list(set(col_cat_list) & set(metadata_vars))
        if len(col_cat_list)> 0:

            C_cat=pd.get_dummies(
                metadata_df[col_cat_list].astype(str),
                drop_first=True)
            C_cat.columns = [c.split("_")[0] for c in C_cat.columns]
            C=C_cat.join(metadata_df[metadata_vars].drop(col_cat_list,axis=1))
        else:
            C=metadata_df[metadata_vars]
        
        #order profiles and C in the same manner, remove samples with missing values
        common_samples=list(set(data_processed_df[idx].dropna().index) & set(C.index))
        '''
        if len(common_samples)<len(data_processed_df.index):
            log_string += "###WARNING###: Not all samples in profiles are found/annotated in metadata. Incomplete samples are being dropped.\n"
        if len(common_samples) < len(C):
            log_string += "###WARNING###: Not all samples in metadata are found in feature data. Incomplete samples are being dropped.\n"
        '''
            
        C=C.loc[common_samples]
        #data_df_sel=data_processed_df.loc[common_samples]

        # C_f -> independent variables for adjustment fit
        if fitsamples_idx_list is not None:
            C_f=C.loc[list(set(common_samples) & set(fitsamples_idx_list))]
        else:
            C_f=C
        #define indices for fit
        sample_index=C_f.index
        
        if len(sample_index) < min_samples:
            log_string += f"\t\t\t\t\t###WARNING###: not enough measured values (less than 10 samples) for given ID. ID {idx} skipped.\n"
            continue
        
        #create a linear model with significantly associated/correlated metadata variables
        
        # z-scale independent variables
        # transformer will be fitted on defined samples only and applied 
        # to all data
        scaler=StandardScaler()
        scaler.fit(C_f)
        C_fz=scaler.transform(C_f)
        C_z=scaler.transform(C)

        if model_string=="ols":
            estimator=linear_model.LinearRegression()
        elif model_string=="ridge":
            estimator=linear_model.RidgeCV(
                alphas=alphas_models)           
        elif model_string=="lasso":
            estimator=linear_model.LassoCV(
                alphas=alphas_models, random_state=random_state)
        elif model_string=="elasticnet":
            estimator=linear_model.ElasticNetCV(
                alphas=alphas_models, random_state=random_state)
        elif model_string=="huber":
            estimator=linear_model.HuberRegressor()
        elif model_string=="bayesian":
            estimator=linear_model.BayesianRidge()
        elif model_string=="logistic":
            estimator=linear_model.LogisticRegression()
        # elif model_string=="GLM":
        #     pass

        median_abundance=data_processed_df.loc[common_samples, idx].median()
        if model_string in ["ols", "ridge", "lasso", "elasticnet", "huber", "bayesian"]:
            estimator.fit(C_fz, data_processed_df.loc[sample_index, idx])
            corr_df.loc[idx,"R2"] = estimator.score(C_fz,data_processed_df.loc[sample_index, idx])
            coefs = estimator.coef_
            n = estimator.intercept_  
            
            # apply correction
            for d in range(len(coefs)):
                data_processed_df.loc[common_samples, idx]=data_processed_df.loc[common_samples, idx] - coefs[d] * C_z[:,d]
            data_processed_df.loc[common_samples, idx] = data_processed_df.loc[common_samples, idx] - n + median_abundance

        elif model_string=="loess":
            estimator=loess.loess(C_fz, data_processed_df.loc[sample_index, idx])
            try:
                estimator.fit()
                pred_values=estimator.predict(C_z).values
                #pred_values_f=estimator.predict(C_fz).values
                corr_df.loc[idx,"R2_"+";".join(metadata_vars)]= metrics.r2_score(C_fz,data_processed_df.loc[sample_index, idx]) 
                data_processed_df.loc[common_samples, idx] = data_processed_df.loc[common_samples, idx] - pred_values + median_abundance
            except ValueError:
                log_string += f"###WARNING###: Fitting did not work (probably too similar outcomes), ID {idx} skipped.\n"
                continue
    
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S") +"\tAdjustment finished!"
    
    return data_processed_df, log_string, corr_df

def adjust_confounders_batchwise(
    data_df,
    metadata_df, 
    model_string, batch_level, random_state, 
    covariates, covariates_dtypes_list,
    alpha_significance=0.05, alphas_models=[-3,3,7],
    fitsamples_idx_list=None, correct_one_by_one=True, 
    correct_over_same_ids_per_batch=False):

    """
    Applies "adjust_confounders" function to batches of data, either 
    for all defined covariates at once or for every covariate
    separately. 

    author: BZ, CD

    Parameters
    ----------
    data_df : pandas.DataFrame
        DataFrame containing the input feature data.
        Rows are samples; columns are features (e.g., gene expression).

    metadata_df : pandas.DataFrame
        DataFrame containing sample-level metadata.
        Must include columns specified in `covariates` and the batch identifier column.

    model_string : str
        The regression model used for confounder adjustment. Options:
        - 'ols': Ordinary Least Squares
        - 'ridge': Ridge Regression (L2 regularization)
        - 'lasso': Lasso Regression (L1 regularization)
        - 'elasticnet': ElasticNet Regression (L1 + L2)
        - 'huber': Huber Regressor (robust to outliers)
        - 'bayesian': Bayesian Ridge Regression
        - 'logistic': Logistic Regression (classification, not recommended for continuous data)
        - 'loess': LOESS (non-parametric, local smoothing)

    batch_level : list of str
        A list containing the column name(s) in `metadata_df` representing batch identifiers.

    random_state : int
        Random seed for reproducibility in models that involve randomness (e.g., LassoCV).

    covariates : list of str
        Metadata column names to be used as confounders for adjustment.

    covariates_dtypes_list : list of str
        List indicating the data type of each covariate. Options:
        - 'nominal'
        - 'ordinal'
        - 'continuous'

    alpha_significance : float, optional (default=0.05)
        Threshold for statistical significance. Covariates with p-values below this 
        threshold are considered for adjustment.

    alphas_models : list-like, optional (default=[-3, 3, 7])
        Range (as powers of 10) for regularization strength in models like Ridge, 
        Lasso, or ElasticNet. Passed to `np.logspace()`.

    fitsamples_idx_list : list-like, optional
        Subset of sample indices to use for model fitting. If None, all samples 
        are used for fitting.

    correct_one_by_one : bool, optional (default=True)
        If True, corrects for each covariate independently. 
        If False, corrects for all covariates simultaneously.
        
    correct_over_same_ids_per_batch : bool, optional (default=True)
        If True, corrects all IDs which are associated with the covariate in any batch.
        If False, correct only IDs which are associated with the covariate per batch.

    Returns
    -------
    data_processed_df : pandas.DataFrame
        The adjusted feature matrix after removing confounding effects.

    log_string : str
        Log of the adjustment process, including messages, warnings, or errors.

    corr_df : pandas.DataFrame
        DataFrame summarizing associations between covariates and features 
        (e.g., R², p-values, correlations) before/after adjustment.
    """

    data_processed_df = data_df.copy()
    log_string = ""
    #select samples used to fit the model
    if fitsamples_idx_list is None:
        fitsamples_idx_list = list(data_df.index)

    corr_dfs=list()
    processed_df_list=list()
    
    #get batches
    batch_levels = sorted(metadata_df[batch_level[0]].unique())
    
    if correct_one_by_one: #correct for one covariate at a time
        for covariate, covariate_dtype in zip(covariates, covariates_dtypes_list):
            processed_batch_df_list = list()
            corr_dfs_covariate_list = list()
            #identify associated ids per batch
            if correct_over_same_ids_per_batch:
                corr_dfs_association = list()
                for batch in batch_levels:
                    #subset fitsamples to fitsamples in given batch
                    idx_batch_list = list(metadata_df[metadata_df[batch_level[0]]==batch].index)
                    fitsamples_idx_batch_list = list(metadata_df.loc[list(set(idx_batch_list)&set(fitsamples_idx_list))].index)
                    corr_df_batch, log_str = get_ids_for_correction(data_df = data_processed_df,
                                                                    metadata_df = metadata_df,
                                                                    adjust_cols_list = [covariate],
                                                                    adjust_cols_dtypes_list = [covariate_dtype],
                                                                    fitsamples_idx_list = fitsamples_idx_batch_list,
                                                                    alpha_significance = alpha_significance)
                    corr_df_batch.columns = [col + "_batch"+str(batch) for col in corr_df_batch.columns]
                    corr_dfs_association.append(corr_df_batch)
                corr_df_covariate = pd.concat(corr_dfs_association,axis=1)
                #take min over corr_dfs to include all precursors significantly associated with covariates at any batch
                corr_df_covariate["pval_" + covariate] = corr_df_covariate[list(filter(lambda x: "pval_" + covariate in x, corr_df_covariate.columns))].min(axis=1)
            else:
                corr_df_covariate = None
            #batchwise correction
            for batch in tqdm(batch_levels):
                idx_batch_list = list(metadata_df[metadata_df[batch_level[0]]==batch].index)
                fitsamples_idx_batch_list = list(metadata_df.loc[list(set(idx_batch_list)&set(fitsamples_idx_list))].index)
                if correct_over_same_ids_per_batch:
                    corr_df_batch = corr_df_covariate[list(filter(lambda x: "_batch"+str(batch) in x, corr_df_covariate.columns)) + ["pval_" + covariate]]
                else:
                    corr_df_batch = None
                data_processed_batch_df, log_string_batch, corr_df_batch =  adjust_confounders(data_df = data_processed_df.loc[idx_batch_list],
                                                                                               metadata_df = metadata_df.loc[idx_batch_list],
                                                                                               model_string = model_string,
                                                                                               adjust_cols_list = [covariate],
                                                                                               adjust_cols_dtypes_list = [covariate_dtype],
                                                                                               random_state = random_state,
                                                                                               fitsamples_idx_list = fitsamples_idx_batch_list,
                                                                                               alpha_significance = alpha_significance,
                                                                                               alphas_models = np.logspace(*alphas_models),
                                                                                               corr_df = corr_df_batch)
                log_string += log_string_batch
                if correct_over_same_ids_per_batch:
                    corr_df_batch.rename(columns={"R2":"R2_batch"+str(batch)},inplace=True)
                    corr_dfs_covariate_list.append(corr_df_batch.drop("pval_" + covariate,axis=1))
                else:
                    corr_df_batch.columns = [col + "_batch"+str(batch) for col in corr_df_batch.columns]
                    corr_dfs_covariate_list.append(corr_df_batch)
                processed_batch_df_list.append(data_processed_batch_df)
            #combine corrected data over batches
            data_processed_df = pd.concat(processed_batch_df_list,axis=0)
            
            corr_df_covariate = pd.concat(corr_dfs_covariate_list,axis=1)
            if correct_over_same_ids_per_batch:
                #take min over corr_dfs to include all precursors significantly associated with covariates at any batch
                corr_df_covariate["pval_" + covariate] = corr_df_covariate[list(filter(lambda x: "pval_" + covariate in x, corr_df_covariate.columns))].min(axis=1)
            #combine_corr_dfs
            corr_dfs.append(corr_df_covariate)
        corr_df = pd.concat(corr_dfs,axis=1)
    else: #correct for all covariates at once
        #identify associated ids per batch
        if correct_over_same_ids_per_batch:
            #create a list with corr_df for every batch
            corr_dfs_association = list()
            for batch in batch_levels:
                #subset fitsamples to fitsamples in given batch
                idx_batch_list = list(metadata_df[metadata_df[batch_level[0]]==batch].index)
                fitsamples_idx_batch_list = list(metadata_df.loc[list(set(idx_batch_list)&set(fitsamples_idx_list))].index)
                corr_df_batch, log_str = get_ids_for_correction(data_df = data_processed_df,
                                                                metadata_df = metadata_df,
                                                                adjust_cols_list = covariates,
                                                                adjust_cols_dtypes_list = covariates_dtypes_list,
                                                                fitsamples_idx_list = fitsamples_idx_batch_list,
                                                                alpha_significance = alpha_significance)
                corr_df_batch.columns = [col + "_batch"+str(batch) for col in corr_df_batch.columns]
                corr_dfs_association.append(corr_df_batch)
            #concatenate all corr_dfs over batches
            corr_pre_df = pd.concat(corr_dfs_association,axis=1)
            #take min over corr_dfs to include all precursors significantly associated with covariates at any batch
            corr_pre_df["pval"] = corr_pre_df[list(filter(lambda x: "pval" in x, corr_pre_df.columns))].min() #HANDLE DIFFERENT COVARIATES #TODO
        else:
            corr_pre_df = None
        #batchwise correction
        for batch in tqdm(batch_levels):
            #subset fitsamples to fitsamples in given batch
            idx_batch_list = list(metadata_df[metadata_df[batch_level[0]]==batch].index)
            fitsamples_idx_batch_list = list(metadata_df.loc[list(set(idx_batch_list)&set(fitsamples_idx_list))].index)
            if len(fitsamples_idx_batch_list)==0:
                warning_batch_level_str = "###ERROR###: No correction variable defined.\n"
                log_string += warning_batch_level_str
                log_string += "###ERROR###: No samples for model-fitting present in batch "+str(batch)+"\n\t\t\t\tPlease adjust fitsamples argument.\n" 
                return data_processed_df, log_string, corr_dfs
            
            #call correction
            (
                data_processed_batch_df, log_string_batch, 
                corr_df_batch
            ) =  adjust_confounders(
                data_df = data_processed_df.loc[idx_batch_list],
                metadata_df = metadata_df.loc[idx_batch_list],
                model_string = model_string,
                adjust_cols_list = covariates,
                adjust_cols_dtypes_list = covariates_dtypes_list,
                random_state = random_state,
                fitsamples_idx_list = fitsamples_idx_batch_list,
                alpha_significance = alpha_significance,
                alphas_models = np.logspace(*alphas_models),
                corr_df = corr_pre_df
                )  
            
            log_string += log_string_batch
            processed_df_list.append(data_processed_batch_df)
            
            #save corr_df_batch to list
            corr_df_batch.columns = [col + "_batch"+str(batch) for col in corr_df_batch.columns]
            corr_dfs.append(corr_df_batch)
    
        # #concatenate all corrected profiles
        data_processed_df=pd.concat(processed_df_list,axis=0)
        
        #concatenate all corr_dfs over batches
        corr_df = pd.concat(corr_dfs,axis=1)
        
        if correct_over_same_ids_per_batch:
            #add r2 values
            corr_df = corr_pre_df.join(corr_df[list(filter(lambda x: "R2" in x, corr_df.columns))]) #CHECK IF THIS IS NECCESSARY AT ALL #TODO
        
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tBatchwise Adjustment finished!\n"
    
    return data_processed_df, log_string, corr_df    

#%%

if __name__ == "__main__":
    
    module = "CORRECTION"
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
    
    (
        data_df, parameters_dict, id_mapping_df, 
        seed, log_string, metadata_df, masked_df, 
        mask_description_df, out_dir_fp, outname, log_file, counter
        ) = helper.init_module(
            module = module, 
            get_metadata=True, 
            get_id_mapping=True, 
            get_mask=True, 
            return_data_fp_only=False,
            skip_parser=skip_parser,
            data_df_fp=data_df_fp,
            param_fp=param_fp,
            out_dir_fp=out_dir_fp,
            meta_df_fp=meta_df_fp,
            log_file=log_file,
            idmap_fp=idmap_fp,
            counter=counter,
            outname=outname)


    # Check if loading was successful, else pass error to log and exit   
    helper.logger(log_file,log_string)
    if metadata_df is None:
        helper.logger(log_file,"###ERROR###: Metadata loading failed.\n")
        if exit_module:
            sys.exit("###ERROR###: Metadata loading failed.")
    if "###WARNING###" in log_string:
        helper.logger(log_file,"###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])
        if exit_module:
            sys.exit("###ERROR###: Metadata does not contain specified metadata_filename_column " + parameters_dict["metadata_filename_column"])

    # reset log string
    log_string = ""
    # set index in metadata
    metadata_df.set_index(parameters_dict["metadata_sample_column"], inplace=True, drop=False)
   
    # -------------------------------------------------------------------------#
    # Check Compatibility of feature data and meta data                        #
    # -------------------------------------------------------------------------#
    
    # check if batchcolumns in metafile
    if parameters_dict["correct_method"] == "harmony" or parameters_dict["correct_method"] == "combat":
        log_string_batchcol = ""

        for batch_col in parameters_dict["batch_cols_meta"]:
            if batch_col not in metadata_df.columns:
                log_string_batchcol += f"###ERROR###: Batch Column {batch_col} not in metafile columns!"
        if "###ERROR" in log_string_batchcol:
            log_string += log_string_batchcol
            helper.logger(log_file, log_string)
            if exit_module:
                sys.exit(1)

    # check if data and meta file consist of same samples
    # only keep intersect on index
    common_idx_list = list(set(data_df.columns) & set(metadata_df.index))  # copy here
    drop_samples_data_lst = list(set(data_df.columns) - set(metadata_df.index))
    drop_samples_meta_lst = list(set(metadata_df.index) - set(data_df.columns))
    
    n_samples_data_int = len(data_df.columns)
    n_samples_meta_int = metadata_df.shape[0]
    data_df = data_df.loc[:, common_idx_list]
    metadata_df = metadata_df.loc[common_idx_list]
    # copy end

    '''
    # if discrepancy between data and metadata samples, log them
    flag_sample_intersect_int = False
    warning_intersect_data_str = ""
    if n_samples_data_int != data_df.shape[1]:
        flag_sample_intersect_int = True
        n_final_samples_data_int = n_samples_data_int - data_df.shape[1]
        warning_intersect_samples_str = \
                f"###WARNING###: {n_final_samples_data_int} sample(s) did not appear in meta data:"\
                + "\n\t\t\t\tFollowing sample(s) has/have been dropped for analysis:\n\t\t\t\t" \
                + "\n\t\t\t\t".join(drop_samples_data_lst) \
                + "\n"
        log_string += warning_intersect_samples_str
    if n_samples_meta_int != metadata_df.shape[0]:
        flag_sample_intersect_int = True
        n_final_samples_meta_int = n_samples_meta_int - metadata_df.shape[0]
        warning_intersect_samples_str = \
                f"###WARNING###: {n_final_samples_meta_int} sample(s) did not appear in precursor data:"\
                + "\n\t\t\t\tFollowing sample(s) has/have been dropped for analysis:\n\t\t\t\t" \
                + "\n\t\t\t\t".join(drop_samples_meta_lst) \
                + "\n"
        log_string += warning_intersect_samples_str
    
    if flag_sample_intersect_int:
        log_string += time.strftime('%Y-%m-%d\t%H:%M:%S') \
                +f"\t{metadata_df.shape[0]} sample(s) remain for analysis.\n"
    '''
    
    # Warn user that some samples do not have a batch annotation
    if parameters_dict["correct_method"] == "harmony" or parameters_dict["correct_method"] == "combat":
        for batch_col in parameters_dict["batch_cols_meta"]:
            samples_na_list = []
            if metadata_df[batch_col].isna().sum() > 0:
                count_na_int = metadata_df[batch_col].isna().sum()
                samples_na_list = list(
                    metadata_df[metadata_df[batch_col].isna()].index)
                warning_na_meta_str = \
                    f"###ERROR###: {count_na_int} sample(s) do(es) not have an {batch_col} value:"\
                    + "\n\t\t\t\t".join(samples_na_list) + "\n"
                if "###ERROR" in log_string_batchcol:
                    log_string += log_string_batchcol
                    helper.logger(log_file, log_string)
                    if exit_module:
                        sys.exit(1)
    
                # drop samples
                if len(samples_na_list) > 0:
                    drop_samples = set(samples_na_list) & set(data_df.index)
                    data_df = data_df.drop(drop_samples, axis=0)
                    metadata_df = metadata_df.drop(drop_samples, axis=0)

        helper.logger(log_file, log_string)

        # Exit module if there are multiple levels in batch_cols of meta file
        for b in parameters_dict["batch_cols_meta"]:
            if len(np.unique(metadata_df[b])) == 1:
                warning_batch_level_str = \
                    f"###ERROR###: Batch Level {b} does only have 1 level:{np.unique(metadata_df[b])[0]} value:"\
                    + "\n\t\t\t\tMultiple levels are required to perform batch correction"  + "\n"
                log_string += warning_batch_level_str
                helper.logger(log_file, log_string)
                if exit_module:
                    sys.exit(1)

    # Add check if some precursors are not present in every batch 
    if parameters_dict["correct_method"] in ["combat", "harmony"]:
        if parameters_dict["check_feature_presence_batch"]:
            batch_presence_th = parameters_dict["presence_batch_threshold"]
            
            flag_presence_batch = False
            log_str_presence = ""

            presence_dfs = dict().fromkeys(
                parameters_dict["batch_cols_meta"], dict())
            missing_precursors_batch = set()
            
            
            # Check if presence is below th in any batch
            for batch in parameters_dict["batch_cols_meta"]:
                presence_batch_df = (
                    data_df.T
                    .notna()
                    .groupby(metadata_df[batch])
                    .mean()
                )
                if (presence_batch_df<batch_presence_th).sum(axis=1).sum() > 0:
                    log_str_presence += f"###ERROR###: Not every precursor has a presence of >{batch_presence_th*100}% in every batch.\nThe following number of precursors are missing per batch:\n"
                    for batch_level in presence_batch_df.index:
                        n_missing_batch_level = (
                            presence_batch_df<batch_presence_th).sum(axis=1).loc[batch_level]
                        log_str_presence += f"\t {batch} {batch_level} :{n_missing_batch_level} \t\t\n"
                    presence_dfs[batch][batch_level] = list(
                        data_df.index[(presence_batch_df<batch_presence_th).loc[batch_level]==True]
                        )
                    missing_precursors_batch = missing_precursors_batch.union(set(presence_dfs[batch][batch_level]))
                    
                    flag_presence_batch = True
            missing_precursors_batch = list(missing_precursors_batch)


            if flag_presence_batch:
                helper.logger(log_file, log_str_presence)
            else:
                helper.logger(log_file, f"All precursor have presence of > {batch_presence_th*100}% in all batches.")



            if flag_presence_batch:
                presence_doc_df = (
                    pd.DataFrame(
                        [
                            (precursor, f"{batch}_{batch_level}")
                            for batch, levels in presence_dfs.items()
                            for batch_level, precursors in levels.items()
                            for precursor in precursors
                        ],
                        columns=["Precursor.Id", "include_correction_batch"],
                    )
                    .groupby("Precursor.Id", as_index=True)
                    .agg({"include_correction_batch":lambda x: ", ".join(sorted(set(x)))})
                )
                
                # manipulate id-mapping by adding missingngess in batch
                id_mapping_df["include_correction"] = True
                id_mapping_df = pd.concat([id_mapping_df, presence_doc_df], axis=1)
                id_mapping_df.loc[
                    missing_precursors_batch, "include_correction"] = False

                id_mapping_df.to_csv(
                    os.path.join(out_dir_fp, "id_mapping.tsv"), index=False, sep="\t")
                sys.exit("###ERROR###: Precursor presence in some batches below threshold!")
    
    # transpose data
    data_df = data_df.T
    # reset list of fit samples
    fitsamples_idx_list = None
    
    # -------------------------------------------------------------------------#
    # Perform Correction                                                 #
    # -------------------------------------------------------------------------#
    
    try:
                if parameters_dict["correct_method"] == "harmony":
            
                    # translate harmony parameter from dict 
                    harmony_kwargs_dict = {
                            "_".join(k.split("_")[1:]):v for k,v in parameters_dict.items()
                            if k.startswith("harmony")}
            
                    data_processed_df, log_string_harmony = batchcorrect_harmony(
                        data_df, metadata_df, 
                        parameters_dict["batch_cols_meta"],
                        seed=seed, 
                        **harmony_kwargs_dict
                    )
                    log_string += log_string_harmony
            
                elif parameters_dict["correct_method"] == "combat":
            
                    # Check if only one batch col is defined!
                    if len(parameters_dict["batch_cols_meta"]) > 1:
                        info_batch_str = \
                                "###ERROR###: COMBAT is only implemented for one batch variable.\n"\
                                +f"\n\t\t\t\t{len(parameters_dict['batch_cols_meta'])} levels are provided: {parameters_dict['batch_cols_meta']}"
                        log_string += info_batch_str
                        helper.logger(log_file,log_string)
                        if exit_module:
                            sys.exit(1)
            
                    # translate combat parameter from dict 
                    combat_kwargs_dict = {
                            "_".join(k.split("_")[1:]):v for k,v in parameters_dict.items()
                            if k.startswith("combat")}
                    if combat_kwargs_dict["fitsamples_meta_column"] is not None and combat_kwargs_dict["fitsamples_meta_level"] is not None:
                        fit_col = combat_kwargs_dict["fitsamples_meta_column"]
                        fit_level = combat_kwargs_dict["fitsamples_meta_level"]
                        
                        fitsamples_idx_list = list(metadata_df[metadata_df[fit_col]==fit_level].index)
                        info_fitsamples_str = \
                                time.strftime("%Y-%m-%d%t%H:%M:%S") \
                                +f"\tCOMBAT will be fitted on {len(fitsamples_idx_list)} samples."\
                                +f"\n\t\t\t\tLevel {fit_level} in meta column {fit_col} will be used to fit COMBAT Model."  + "\n"
                        log_string += info_fitsamples_str
            
                        # Check whether the same batch levels are defined in samples
                        # fitting combat and the ones and the rest of the dataset
            
                        levels_fitsamples = set(
                            metadata_df.loc[
                                fitsamples_idx_list, parameters_dict["batch_cols_meta"][0]])
                        levels_data = set(
                            metadata_df[parameters_dict["batch_cols_meta"][0]])
                        if len(levels_data) != len(levels_fitsamples):
                            info_fitsamples_str = \
                                f"###WARNING###: Uneven numbers of batch levels in whole dataset (n={len(levels_data)}) and samples selected to fit COMBAT model (n={len(levels_fitsamples)})."\
                                +f"\n\t\t\t\tLevel(s) {list(levels_data-levels_fitsamples)} do not appear in samples to fit model and samples will be discarded"  + "\n"
                            log_string += info_fitsamples_str
                        # Remove samples of batch levels do not appear in fit samples
                        metadata_df = metadata_df[
                            ~metadata_df[parameters_dict["batch_cols_meta"][0]].isin(
                                list(levels_data-levels_fitsamples))]
                        data_df = data_df.loc[metadata_df.index]


                    if parameters_dict["combat_cols_adjust_dtypes"] is None:
                                parameters_dict["combat_cols_adjust_dtypes"] = []
            
                    out = batchcorrect_combat(
                        data_df = data_df,# possibly data_imputed_df
                        metadata_df = metadata_df,
                        batch_key = parameters_dict["batch_cols_meta"][0],
                        fitsamples_idx_list=fitsamples_idx_list,
                        cov_preserve_list = parameters_dict["combat_cols_preserve_meta"],
                        cov_adjust_list = parameters_dict["combat_cols_adjust_meta"],
                        cov_adjust_list_dtypes = parameters_dict["combat_cols_adjust_dtypes"]
                    )
                    if isinstance(out, str) and "ERROR" in out:
                        log_string_combat = out 
                        log_string += log_string_combat
                        helper.logger(log_file, log_string)
                        if exit_module:
                            sys.exit(1)
                    else:
                        (data_processed_df, log_string_combat) = out
                        log_string += log_string_combat
            
                elif parameters_dict["correct_method"] == "adjust_confounders":
            
                    # check if number of cols and datatypes have same length
                    if not len(parameters_dict["adjust_cols"]) == len(parameters_dict["adjust_cols_dtypes"]):
                        log_string += "###ERROR###: Number of columns for adjustment and number of datatypes are different.\n"
                        helper.logger(log_file, log_string)
                        if exit_module:
                            sys.exit(1)
            
                    
                    valid_dtypes_list = ["continuous", "nominal", "ordinal"]
                    invalid_dtypes_list = []
                    for dt in parameters_dict["adjust_cols_dtypes"]:
                        if dt not in valid_dtypes_list:
                            invalid_dtypes_list.append(dt)
                    
                    if len(invalid_dtypes_list) > 0:
                        log_string += f"###ERROR###: The following datatypes are not valid: {invalid_dtypes_list}.\n"
                        helper.logger(log_file, log_string)
                        if exit_module:
                            sys.exit(1)
            
            
                    # Check for missing values in columns for adjustment
                    for var in parameters_dict["adjust_cols"]:
                        if metadata_df[var].isna().any():
                            n_na = metadata_df[metadata_df[var].isna()].shape[0]
                            log_string += f"###ERROR###: Column {var} contains n={n_na} missing values.\n"
                    if "###ERROR###" in log_string:
                        helper.logger(log_file, log_string)
                        if exit_module:
                            sys.exit(1)
            
            
            
                    # Check if datatype can be casted to float or string
                    for var, dt in zip (
                        parameters_dict["adjust_cols"], 
                        parameters_dict["adjust_cols_dtypes"]):
                        if dt == "continuous":
                            try:
                                metadata_df[var] = metadata_df[var].astype(float)
                            except:
                                log_string += f"###ERROR###: Column {var} could not be casted to {dt}.\n"
            
                        elif dt == "ordinal":
                            try:
                                metadata_df[var] = metadata_df[var].astype(int)
                            except:
                                log_string += f"###ERROR###: Column {var} could not be casted to {dt}.\n"
                        
                        elif dt == "nominal":
                            #nominal variable should be casted to string -> object in pandas
                            try:
                                metadata_df[var] = metadata_df[var].astype(str)
                            except:
                                log_string += f"###ERROR###: Column {var} could not be casted to {dt}.\n"
                        else:
                            log_string += "###ERROR###: Datatype not valid.\n"                
            
                    # Check if requested model is implemented
                    if parameters_dict["adjust_model"].lower() not in [
                            "ols","ridge","lasso","elasticnet","huber","bayesian","loess"]:
                            
                            log_string += f"###ERROR###: Model {parameters_dict['adjust_model']} not implemented, no correction applied.\n"                
                    
                    #if adjustment model is LOESS, switch adjust_covariates_onebyone on because LOESS correction cannot deal with several variables
                    if parameters_dict["adjust_model"] == "loess" and not parameters_dict["adjust_covariates_onebyone"]:
                        parameters_dict["adjust_covariates_onebyone"] = 1
                        log_string += "###WARNING###: Variable adjust_covariates_onebyone set to 1 for adjust_model LOESS, LOESS cannot be applied to correct for several variables at once.\n"
                    
                    if parameters_dict["adjust_fitsamples_meta_column"] is not None and parameters_dict["adjust_fitsamples_meta_level"] is not None:
                        
                        fit_col = parameters_dict["adjust_fitsamples_meta_column"]
                        fit_level = parameters_dict["adjust_fitsamples_meta_level"]
                        fitsamples_idx_list = list(metadata_df[metadata_df[fit_col]==fit_level].index)
                        info_fitsamples_str = \
                                time.strftime("%Y-%m-%d%t%H:%M:%S") \
                                +f"\tAdjustment models will be fitted on {len(fitsamples_idx_list)} samples."\
                                +f"\n\t\t\t\tLevel {fit_level} in meta column {fit_col} will be used to adjust samples."  + "\n"
                        log_string += info_fitsamples_str
            
                        # Check whether the same levels are defined in samples #DO WE NEEDED? Does not work for continuous variables anyway #???
                        """
                        levels_fitsamples = set(
                            metadata_df.loc[
                                fitsamples_idx_list, parameters_dict["adjust_cols"][0]])
                        levels_data = set(
                            metadata_df[parameters_dict["adjust_cols"][0]])
                        if len(levels_data) != len(levels_fitsamples):
                            info_fitsamples_str = \
                                time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                                +f"###WARNING###: Uneven numbers of batch levels in whole dataset (n={len(levels_data)}) and samples selected to fit adjustment model (n={len(levels_fitsamples)}).\n"\
                                +f"Level(s) {list(levels_data-levels_fitsamples)} do not appear in samples to fit model and samples will be discarded"  + "\n"
                            log_string += info_fitsamples_str
                        # Remove samples of batch levels do not appear in fit samples
                        metadata_df = metadata_df[
                            ~metadata_df[parameters_dict["adjust_cols"][0]].isin(
                                list(levels_data-levels_fitsamples))]
                        data_df = data_df.loc[metadata_df.index]
                        """
                        
                    if "ERROR" in log_string:
                        # Exit module because of error
                        helper.logger(log_file,log_string)
                        if exit_module:
                            sys.exit(1)
                    
                    else:
                        if parameters_dict["adjust_covariates_onebyone"]:   
                            data_processed_df = data_df
                            corr_df, log_str = get_ids_for_correction(data_df=data_processed_df,
                                                                      metadata_df=metadata_df,
                                                                      adjust_cols_list=parameters_dict["adjust_cols"],
                                                                      adjust_cols_dtypes_list=parameters_dict["adjust_cols_dtypes"],
                                                                      fitsamples_idx_list=fitsamples_idx_list,
                                                                      alpha_significance=parameters_dict["adjust_alpha"])
                            log_string += log_str
                            for c in range(len(parameters_dict["adjust_cols"])):
                                data_processed_df, log_string_adjust, corr_df = adjust_confounders(
                                    data_df = data_processed_df,
                                    metadata_df = metadata_df,
                                    model_string = parameters_dict['adjust_model'],
                                    adjust_cols_list = [parameters_dict["adjust_cols"][c]],
                                    adjust_cols_dtypes_list = [parameters_dict["adjust_cols_dtypes"][c]],
                                    random_state = seed,
                                    fitsamples_idx_list = fitsamples_idx_list,
                                    alpha_significance = parameters_dict["adjust_alpha"],
                                    alphas_models = np.logspace(*parameters_dict["adjust_alphas_models"]),
                                    corr_df = corr_df)
                                log_string += log_string_adjust
                                
                        else:
                            data_processed_df, log_string_adjust, corr_df = adjust_confounders(
                                data_df = data_df,
                                metadata_df = metadata_df,
                                model_string = parameters_dict['adjust_model'],
                                adjust_cols_list = parameters_dict["adjust_cols"],
                                adjust_cols_dtypes_list = parameters_dict["adjust_cols_dtypes"],
                                random_state = seed,
                                fitsamples_idx_list = fitsamples_idx_list,
                                alpha_significance = parameters_dict["adjust_alpha"],
                                alphas_models = np.logspace(*parameters_dict["adjust_alphas_models"]),
                                corr_df = None)
                        
                        #save dataframe with association pvalues, effect sizes and R2 after correction
                        corr_df.to_csv(os.path.join(out_dir_fp,"precursor_associations.tsv"),sep="\t")
            
                elif parameters_dict["correct_method"] == "adjust_confounders_batchwise":
                    
                    if parameters_dict["adjust_cols"] is None:
                        warning_batch_level_str = \
                            "###ERROR###: No correction variable defined \n"
                        log_string += warning_batch_level_str
                        helper.logger(log_file, log_string)
                        if exit_module:
                            sys.exit(1)
                    
                    #if adjustment model is LOESS, switch adjust_covariates_onebyone on because LOESS correction cannot deal with several variables
                    if parameters_dict["adjust_model"] == "loess" and not parameters_dict["adjust_covariates_onebyone"]:
                        parameters_dict["adjust_covariates_onebyone"] = 1
                        log_string += "###WARNING###: Variable adjust_covariates_onebyone set to 1 for adjust_model LOESS, LOESS cannot be applied to correct for several variables at once.\n"
            
            
                    fitsamples_idx_list = None
                    if parameters_dict["adjust_fitsamples_meta_column"] is not None and parameters_dict["adjust_fitsamples_meta_level"] is not None:
                        
                        fit_col = parameters_dict["adjust_fitsamples_meta_column"]
                        fit_level = parameters_dict["adjust_fitsamples_meta_level"]
                        #check if fit_col is in metadata_df
                        if fit_col not in metadata_df.columns:
                            log_string += \
                            f"###ERROR###: Column '{fit_col}' as defined by 'adjust_fitsamples_meta_column' not found in metadata.\n"
                            helper.logger(log_file,log_string)
                            if exit_module:
                                sys.exit(1)
                        # check if fit_value is in metadata_df[fit_col]
                        if fit_level not in metadata_df[fit_col].tolist():
                            log_string += \
                            f"###ERROR###: Value '{fit_level}' as defined by 'adjust_fitsamples_meta_level' not found in column '{fit_col}' in metadata.\n"
                            helper.logger(log_file,log_string)
                            if exit_module:
                                sys.exit(1)
                        fitsamples_idx_list = list(metadata_df[metadata_df[fit_col]==fit_level].index)
                        info_fitsamples_str = \
                                time.strftime("%Y-%m-%d%t%H:%M:%S") \
                                +f"\tCorrection model will be fitted on {len(fitsamples_idx_list)} samples."\
                                +f"\n\t\t\t\tLevel {fit_level} in meta column {fit_col} will be used to fit {parameters_dict['adjust_model']} Model for covariate correction."  + "\n"
                        log_string += info_fitsamples_str
                    
                    data_processed_df, log_string, corr_df = adjust_confounders_batchwise(
                        data_df=data_df,
                        metadata_df=metadata_df, 
                        model_string=parameters_dict["adjust_model"], 
                        batch_level=parameters_dict["batch_cols_meta"], 
                        random_state=seed, 
                        covariates=parameters_dict["adjust_cols"], 
                        covariates_dtypes_list=parameters_dict["adjust_cols_dtypes"],
                        alpha_significance=parameters_dict["adjust_alpha"], 
                        alphas_models=parameters_dict["adjust_alphas_models"],
                        fitsamples_idx_list=fitsamples_idx_list, 
                        correct_one_by_one=parameters_dict["adjust_covariates_onebyone"],
                        correct_over_same_ids_per_batch=parameters_dict["adjust_covariates_over_same_ids_per_batch"])
                    
                    #save dataframe with association pvalues, effect sizes and R2 after correction
                    corr_df.to_csv(os.path.join(out_dir_fp,"precursor_associations.tsv"),sep="\t")
                    
                else:
                    data_processed_df = data_df
                    log_string += f"###WARNING###: Unknown correction method: {parameters_dict['correct_method']} — skipping.\n"
    except Exception as e:
                e_type, e_obj, e_tb = sys.exc_info()
                script_name = os.path.split(e_tb.tb_frame.f_code.co_filename)[1]
                helper.logger(log_file, "###ERROR### Module " +script_name + " ran into an error at line : " + str(e_tb.tb_lineno) +"\n\t\t\t\t" + str(e) + "\n")
                if exit_module:
                    sys.exit("Correction exited with an error.")

    data_processed_df = data_processed_df.T

    helper.close_module(data_processed_df=data_processed_df, 
                        masked_df=masked_df, 
                        mask_description_df=mask_description_df, 
                        log_file=log_file, module_output=log_string,
                        out_dir_fp=out_dir_fp, outname=outname, mask=True, exit_module=exit_module)


