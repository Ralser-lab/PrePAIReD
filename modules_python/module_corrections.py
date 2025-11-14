"""
Input: wide-format data
Output: Corrected wide-format data


- batchcorrection_combat: applies combat batch correction
- batchcorrection_harmony: applies harmony batch correction
- correct_for_confounders: drift correction per id for a given set of confounders which are statistically significant correlated to the id
- adjust_confounders_batchwise: applies drift correction per batch for a set of covariates

General comments:
- No duplicated adjustment in Combat and via linear modeling

# TBD save of corr_dfs in adjust_confounders and adjust_confounders_batchwise missing
# TBD Detailed Module Overview
# Implemented for features as columns and samples as rows

"""
# %%

#general
import os
import sys
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import time

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
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: "+"Start HARMONY batch correction\n")

    data_arr = data_df.to_numpy()
    # Apply harmony algorithm
    data_arr = data_df.to_numpy()
    data_processed_arr = harmonize(
        data_arr, metadata_df, batch_key, random_state=seed,
        **kwargs)
    data_processed_df = pd.DataFrame(
        data_processed_arr, 
        index=data_df.index, columns=data_df.columns)

    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: "\
                      +"Finished HARMONY batch correction")

    return data_processed_df, log_string


def batchcorrect_combat(
    data_df:pd.DataFrame, 
    metadata_df:pd.DataFrame, 
    batch_key:str,
    fitsamples_idx_list=list,
    cov_preserve_list=list,
    cov_adjust_list=list,
    ):
        
        """
        Apply COMBAT batch correction by using the pyComBat package.
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
                "%Y-%m-%d%t%H:%M:%S: "+"Start COMBAT batch correction.\n")

        # effects to preserve during batch correction
        if cov_preserve_list != None:
            x_df = metadata_df[cov_preserve_list]

        # covariates to adjust for during batch correction
        if cov_adjust_list != None:
            c_df = metadata_df[cov_adjust_list]

        # Init combat Model
        combat_model=pycombat.Combat()
        # fit combat model on selected samples
        if fitsamples_idx_list != None:
            
            log_string += time.strftime(
                    "%Y-%m-%d%t%H:%M:%S: "+"Fit COMBAT Class on defined samples.\n")

            Y_fit = data_df.loc[fitsamples_idx_list].to_numpy()
            b_fit = metadata_df[batch_key].loc[fitsamples_idx_list].astype(str).to_numpy()
            X_fit = x_df.loc[fitsamples_idx_list]
            C_fit = c_df.loc[fitsamples_idx_list]

            c_fit = combat_model.fit(
                Y=Y_fit,
                b=b_fit,
                X=X_fit,
                C=C_fit
            )

        else:
            
            log_string += time.strftime(
                    "%Y-%m-%d%t%H:%M:%S: "+"Fit COMBAT Class on all samples.\n")

            c_fit = combat_model.fit(
                Y=data_df.T.to_numpy(),
                b=metadata_df[batch_key].to_numpy(),
                X=x_df, 
                C=c_df
            )


        log_string += time.strftime(
                    "%Y-%m-%d%t%H:%M:%S: "+"Transform data.\n")

        samples_all = data_df.columns
        Y_all = data_df.T.to_numpy()
        b_all = metadata_df[batch_key].loc[samples_all].astype(str).to_numpy()
        X_all = x_df.loc[samples_all]
        C_all = c_df.loc[samples_all]

        assert Y_all.shape[0] == X_all.shape[0] == C_all.shape[0] == b_all.shape[0]

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

        # Check if output df consists of NA's exclusively
        # if so return error in Log and don't return processed df
        if data_processed_df.isna().all().all():
            log_string += time.strftime(
                    "%Y-%m-%d%t%H:%M:%S: "+"###ERROR###: All Values in DF were transformed to NA.\n"\
                    +f"Please filter input dataframe for features with high missing values or low variance."  + "\n")
            return log_string
        else:
            log_string += time.strftime(
                "%Y-%m-%d%t%H:%M:%S: "\
                +"Finished COMBAT batch correction successfully.\n")

            return data_processed_df, log_string
        

def adjust_confounders(
                data_df, metadata_df, model_string, adjust_cols_list,
                adjust_cols_dtypes_list, random_state, fitsamples_idx_list=None,
                alpha_significance=0.05, alphas_models=[-3,3,7]
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
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") + "Start Adjustment.\n"
    
    # use selecterd samples to fit model
    samples_idx_list = data_processed_df.index
    
    if fitsamples_idx_list != None:
        samples_idx_list = fitsamples_idx_list
    
    # collect associations with meta data variables
    corr_df=pd.DataFrame(index=data_processed_df.columns,
                        columns=["R2"])
    corr_df.sort_index(inplace=True)
    # iterate over all features

    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") + f"{len(fitsamples_idx_list)} samples were selected to fit model.\n"

    # iterate over features
    for feat in data_processed_df.columns:
        # iterate over defined columns and assess statistical significance
        metadata_vars = [] # collect variables that show significant associations with feature
        for metadata_var, dt in zip(
            adjust_cols_list,
            adjust_cols_dtypes_list):
            # remove missing values and define commindex in meta and feature data
            idx_meta = metadata_df.loc[samples_idx_list, metadata_var].dropna()
            idx_data = data_processed_df.loc[samples_idx_list, feat].dropna()
            common_idx = list(set(idx_meta.index) & set(idx_data.index))

            samples_idx_test_list= [
                list(metadata_df[metadata_df[metadata_var]==x].index) 
                for x in metadata_df[metadata_var].unique()]
            samples_idx_test_list = [list(set(x) & set(common_idx)) for x in samples_idx_test_list]
            # Check if at least 5 samples are defined per group

            # correlation coefficient to be calculated depending on datatype
            coef_type_str = ""
            if dt == "nominal":
                coef_type_str = "pointbiserialr_"
            elif dt == "ordinal":
                coef_type_str = "spearmanr_"
            elif dt == "continuous":
                coef_type_str = "pearsonr_"
            else:
                pass

            # check how many samples are available per feature
            # do not check association if <1ß
            if len(common_idx) < 5:
                log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    + f"###Warning###: One of the groups in col {metadata_var} consists of less than 5 samples.\n" \
                        f"No Correction will be applied on {feat}.\n"
                corr_df.loc[feat, "pval_"+ metadata_var] = np.nan
                corr_df.loc[feat, coef_type_str+ metadata_var] = np.nan

            # check if ordinal or nomonal variable groups have at 
            # do not check association if <5 samples per group
            elif (dt == "ordinal") or (dt=="nominal"):
                if any([len(x) for x in samples_idx_test_list]) < 5:
                    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    + f"###Warning###: One of the groups in col {metadata_var} consists of less than 5 samples.\n" \
                        f"No Correction will be applied on {feat}.\n"

                    
                    # save nan if sample size was not sufficient
                    corr_df.loc[feat, "pval_"+ metadata_var] = np.nan
                    corr_df.loc[feat, coef_type_str+ metadata_var] = np.nan
            

            else:
                samples_test_list = [data_processed_df.loc[x, feat] for x in samples_idx_test_list]
            
                #test for significant associations/correlations
                if dt=="nominal":
                    if len(metadata_df.loc[common_idx,metadata_var].unique())==2:
                        # cast metadata variable to boolean for pointbiserial
                        if metadata_df.loc[common_idx,metadata_var].dtype not in [bool, float]:
                            levels_var = metadata_df.loc[common_idx,metadata_var].unique()
                            metadata_df[metadata_var] = metadata_df[metadata_var].replace(
                                {levels_var[0]:True, levels_var[1]:False})
                        # Point biserial correlation + one way ANOVA
                        
                        association, pval_corr = stats.pointbiserialr(
                            metadata_df.loc[common_idx,metadata_var],
                            data_processed_df.loc[common_idx, feat])
                        pval_test = stats.f_oneway(*samples_test_list).pvalue
                        pval = pval_test
                        if (pval_test<=alpha_significance and pval_corr>alpha_significance) or (pval_test>alpha_significance and pval_corr<=alpha_significance):
                            log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                            + f"WARNING: p-value from correlation and ANOVA disagree on neglecting null hypothesis for {feat}"
                        corr_df.loc[feat,"pointbiserialr_" + metadata_var]=association
                    elif len(metadata_df.loc[common_idx,metadata_var].unique())<2:
                        #no association possible with one value per variable
                        association=0
                        pval=1
                    else:
                        pval=stats.f_oneway(*samples_test_list).pvalue
                
                elif dt == "ordinal":
                    #Spearman correlation + one way ANOVA
                    association, pval_corr = stats.spearmanr(
                        metadata_df.loc[common_idx, metadata_var],
                        data_processed_df.loc[common_idx, feat])
                    pval_test=stats.f_oneway(*samples_test_list).pvalue
                    pval = pval_test
                    if (pval_test<=alpha_significance and pval_corr>alpha_significance) or (pval_test>alpha_significance and pval_corr<=alpha_significance):
                        log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                        +f"WARNING: p-value from correlation and ANOVA disagree on neglecting null hypothesis for {feat}"
                    corr_df.loc[feat,"spearmanr_" + metadata_var]=association

                elif dt =="continuous":

                    
                    #Pearson correlation
                    association, pval=stats.pearsonr(
                        metadata_df.loc[common_idx, metadata_var],
                        data_processed_df.loc[common_idx, feat])
                    corr_df.loc[feat,"pearsonr_" + metadata_var]=association
                else:
                    pass

                # if p-value of association test below significance value 
                # variable will be stored for later adjustment
                if pval < alpha_significance:
                    metadata_vars.append(metadata_var)
                corr_df.loc[feat, "pval_"+ metadata_var] = pval

            
            if len(metadata_vars)==0:
                continue

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
            
            #order profiles and C in the same manner
            common_samples=list(set(data_processed_df.index) & set(C.index))
    
            if len(common_samples)<len(data_processed_df.index):
                log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                + "###WARNING###: Not all samples in profiles are found/annotated in metadata. Incomplete samples are being dropped."
            if len(common_samples) < len(C):
                log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                +"###WARNING###: Not all samples in metadata are found in feature data. Incomplete samples are being dropped."
            C=C.loc[common_samples]
            data_df_sel=data_processed_df.loc[common_samples]

            # C_f -> independent variables for adjustment fit
            if fitsamples_idx_list is not None:
                C_f=C.loc[common_idx]
            else:
                C_f=C
            #define indices for fit
            sample_index=C_f.index
            
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

            median_abundance=data_processed_df.loc[:, feat].median()
            if model_string in ["ols", "ridge", "lasso", "elasticnet", "huber", "bayesian"]:
                estimator.fit(C_fz, data_processed_df.loc[sample_index, feat])
                corr_df.loc[feat,"R2"] = estimator.score(C_fz,data_processed_df.loc[sample_index, feat])
                coefs = estimator.coef_
                n = estimator.intercept_  
                
                # apply correction
                if fitsamples_idx_list is not None:
                    for d in range(len(coefs)):
                        data_processed_df.loc[:, feat]=data_processed_df.loc[
                            :, feat]- coefs[d] * C_z[:,d]
                    data_processed_df.loc[common_idx, feat] = data_processed_df.loc[
                        common_idx, feat]-n + median_abundance
                else:
                    for d in range(len(coefs)):
                        data_processed_df.loc[:, feat]=data_processed_df.loc[:, feat]-coefs[d]*C_z[:,d]
                    data_processed_df.loc[:, feat]=data_processed_df.loc[:, feat]-n + median_abundance

            elif model_string=="loess":
                estimator=loess.loess(
                    C_fz, data_processed_df.loc[sample_index, feat])
                try:
                    estimator.fit()
                except ValueError:
                    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###WARNING###: Fitting did not work (probably too similar outcomes), feature {feat} skipped."
                    continue
                pred_values=estimator.predict(C_fz).values
                if fitsamples_idx_list is not None:
                    corr_df.loc[feat,"R2"]= metrics.r2_score(
                        data_processed_df.loc[common_idx, feat],pred_values)
                    data_processed_df.loc[common_idx, feat] = data_processed_df.loc[
                        common_idx, feat]-pred_values+median_abundance
                else:
                    corr_df.loc[feat,"R2"] = metrics.r2_score(data_processed_df[feat],pred_values)
                    data_processed_df.loc[:, feat]=data_processed_df.loc[:, feat]-pred_values + median_abundance
    
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") +f"Adjustment finished!"
    
    return data_processed_df, log_string, corr_df


def adjust_confounders_batchwise(
    data_df,
    metadata_df, 
    model_string, batch_level, random_state, 
    covariates, covariates_dtypes_list,
    alpha_significance=0.05, alphas_models=[-3,3,7],
    fitsamples_idx_list=None, correct_one_by_one=True):

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

    if fitsamples_idx_list is None:
        fitsamples_idx_list = list(data_df.index)

    corr_dfs=list()
    batch_levels = sorted(metadata_df[batch_level[0]].unique())

    if correct_one_by_one:
        for covariate, covariate_dtype in zip(
            covariates, covariates_dtypes_list):
            processed_df_list=list()
            
            for batch in batch_levels:
                idx_batch_list = list(
                    metadata_df[
                        metadata_df[batch_level[0]]==batch].index)
                fitsamples_idx_batch_list = list(metadata_df.loc[
                    list(set(idx_batch_list)&set(fitsamples_idx_list))].index)
                (
                    data_processed_batch_df, log_string_batch, 
                    corr_df
                ) =  adjust_confounders(
                    data_df = data_processed_df.loc[idx_batch_list],
                    metadata_df = metadata_df.loc[idx_batch_list],
                    model_string = model_string,
                    adjust_cols_list = [covariate],
                    adjust_cols_dtypes_list = [covariate_dtype],
                    random_state = random_state,
                    fitsamples_idx_list = fitsamples_idx_batch_list,
                    alpha_significance = alpha_significance,
                    alphas_models = np.logspace(*alphas_models)
                    )
                log_string += log_string_batch
                corr_df.columns = [x + " batch"+str(batch) for x in corr_df.columns]
                corr_dfs.append(corr_df)
                processed_df_list.append(data_processed_batch_df)
    else: #correct for all covariates at once
        processed_df_list=list()
        for batch in batch_levels:
            idx_batch_list = list(
                    metadata_df[
                        metadata_df[batch_level[0]]==batch].index)
            fitsamples_idx_batch_list = list(metadata_df.loc[
                list(set(idx_batch_list)&set(fitsamples_idx_list))].index)
            if len(fitsamples_idx_batch_list)==0:
                warning_batch_level_str = \
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###ERROR###: No correction variable defined.\n"
                log_string += warning_batch_level_str
                log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                + "###ERROR###:No samples for model-fitting present in batch "+str(batch)+"\nPlease adjust fitsamples argument.\n" 
                return data_processed_df, log_string, corr_dfs
            (
                data_processed_batch_df, log_string_batch, 
                corr_df
            ) =  adjust_confounders(
                data_df = data_processed_df.loc[idx_batch_list],
                metadata_df = metadata_df.loc[idx_batch_list],
                model_string = model_string,
                adjust_cols_list = covariates,
                adjust_cols_dtypes_list = covariates_dtypes_list,
                random_state = random_state,
                fitsamples_idx_list = fitsamples_idx_batch_list,
                alpha_significance = alpha_significance,
                alphas_models = np.logspace(*alphas_models)
                )  
            log_string += log_string_batch      
            corr_df.columns = [x + " batch"+str(batch) for x in corr_df.columns]
        corr_dfs.append(corr_df)
        processed_df_list.append(data_processed_batch_df)
    
    # #concatenate all corrected profiles
    data_processed_df=pd.concat(processed_df_list)
    corr_dfs=pd.concat(corr_dfs,axis=1)
        
    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
    +f"Batchwise Adjustment finished!"
    return data_processed_df, log_string, corr_dfs    


if __name__ == "__main__":
    
    (
        data_df, parameters_dict, id_mapping_df, 
        seed, module_output, metadata_df, masked_df, 
        mask_description_df, out_dir_fp, outname, log_file, counter
        ) = helper.init_module(
            module = "Correction", 
            get_metadata=True, 
            get_id_mapping=True, 
            get_mask=True, 
            return_data_fp_only=False)

    log_string = time.strftime(
        "%Y-%m-%d%t%H:%M:%S: "+ f"Start batch correction {counter}!\n")

    # -------------------------------------------------------------------------#
    # Check Compatibility of feature and meta data                             #
    # -------------------------------------------------------------------------#
    
    # check if batchcolumns in metafile
    if parameters_dict["correct_method"] == "harmony" or parameters_dict["correct_method"] == "combat":
        for batch_col in parameters_dict["batch_cols_meta"]:
            if batch_col not in metadata_df.columns:
                log_string += f"###ERROR###: Batch Column {batch_col} not in metafile columns!"
                helper.logger(log_file, log_string)
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
         
    flag_sample_intersect_int = False
    warning_intersect_data_str = ""
    if n_samples_data_int != data_df.shape[1]:
        flag_sample_intersect_int = True
        n_final_samples_data_int = n_samples_data_int - data_df.shape[1]
        warning_intersect_samples_str = \
                time.strftime('%Y-%m-%d\t%H:%M:%S: ') \
                +f"###WARNING###: {n_final_samples_data_int} sample(s) did not appear in meta data:"\
                + "\nFollowing sample(s) has/have been dropped for analysis:" \
                + "\n\t"+'\n\t'.join(drop_samples_data_lst) \
                + "\n"
        log_string += warning_intersect_samples_str
    if n_samples_meta_int != metadata_df.shape[0]:
        flag_sample_intersect_int = True
        n_final_samples_meta_int = n_samples_meta_int - metadata_df.shape[0]
        warning_intersect_samples_str = \
                time.strftime('%Y-%m-%d\t%H:%M:%S: ') \
                +f"###WARNING###: {n_final_samples_meta_int} sample(s) did not appear in precursor data:\n"\
                + "\n\t\t Following sample(s) has/have been dropped for analysis:\n" \
                + "\t\t" + "\n\t\t".join(drop_samples_meta_lst) \
                + "\n"
        log_string += warning_intersect_samples_str
    
    if flag_sample_intersect_int:
        log_string += time.strftime('%Y-%m-%d\t%H:%M:%S ') \
                +f"{metadata_df.shape[0]} remain for analysis.\n"
        
    # Warn user that some samples do not have a batch annotation
    if parameters_dict["correct_method"] == "harmony" or parameters_dict["correct_method"] == "combat":
        for batch_col in parameters_dict["batch_cols_meta"]:
            samples_na_list = []
            if metadata_df[batch_col].isna().sum() > 0:
                count_na_int = metadata_df[batch_col].isna().sum()
                samples_na_list = list(
                    metadata_df[metadata_df[batch_col].isna()].index)
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
                    metadata_df = metadata_df.drop(drop_samples, axis=0)

        # Exit module if there are multiple levels in batch_cols of meta file
        for b in parameters_dict["batch_cols_meta"]:
            if len(np.unique(metadata_df[b])) == 1:
                warning_batch_level_str = \
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###ERROR###: Batch Level {b} does only have 1 level:{np.unique(metadata_df[b])[0]} value:"\
                    + "\t\tMultiple levels are required to perform batch correction"  + "\n"
                log_string += warning_batch_level_str
                helper.logger(log_file, log_string)
                sys.exit(1)
    
    
    # -------------------------------------------------------------------------#
    # Perform Batch Correction                                                 #
    # -------------------------------------------------------------------------#
    data_df = data_df.T

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
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###ERROR###: COMBAT is only implemented for one batch variable.\n"\
                    +f"{len(parameters_dict['batch_cols_meta'])} levels are provided: {parameters_dict['batch_cols_meta']}"  + "\n"
            log_string += info_batch_str
            helper.logger(log_file,log_string)
            sys.exit(1)

        # translate combat parameter from dict 
        combat_kwargs_dict = {
                "_".join(k.split("_")[1:]):v for k,v in parameters_dict.items()
                if k.startswith("combat")}
        fitsamples_idx_list = None
        if len(combat_kwargs_dict["fitsamples_meta"]) == 2:
            fit_col, fit_level = combat_kwargs_dict["fitsamples_meta"]
            fitsamples_idx_list = list(metadata_df[metadata_df[fit_col]==fit_level].index)
            info_fitsamples_str = \
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"COMBAT will be fitted on {len(fitsamples_idx_list)} samples.\n"\
                    +f"Level {fit_level} in meta column {fit_col} will be used to fit COMBAT Model."  + "\n"
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
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###WARNING###: Uneven numbers of batch levels in whole dataset (n={len(levels_data)}) and samples selected to fit COMBAT model (n={len(levels_fitsamples)}).\n"\
                    +f"Level(s) {list(levels_data-levels_fitsamples)} do not appear in samples to fit model and samples will be discarded"  + "\n"
                log_string += info_fitsamples_str
            # Remove samples of batch levels do not appear in fit samples
            metadata_df = metadata_df[
                ~metadata_df[parameters_dict["batch_cols_meta"][0]].isin(
                    list(levels_data-levels_fitsamples))]
            data_df = data_df.loc[:, metadata_df.index]


        out = batchcorrect_combat(
            data_df = data_df,# possibly data_imputed_df
            metadata_df = metadata_df,
            batch_key = parameters_dict["batch_cols_meta"][0],
            fitsamples_idx_list=fitsamples_idx_list,
            cov_preserve_list = parameters_dict["combat_cols_preserve_meta"],
            cov_adjust_list = parameters_dict["combat_cols_adjust_meta"],
        )
        if isinstance(out, str) and "ERROR" in out:
            log_string_combat = out 
            log_string += log_string_combat
            sys.exit(1)
        else:
            (data_processed_df, log_string_combat) = out
            log_string += log_string_combat

    elif parameters_dict["correct_method"] == "adjust_confounders":

        # check if number of cols and datatypes have same length
        if not len(parameters_dict["adjust_cols"]) == len(parameters_dict["adjust_cols_dtypes"]):
            log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###ERROR###: Number of columns for adjustment and number of datatypes are different.\n"
            helper.logger(log_file, log_string)
            sys.exit(1)

        
        valid_dtypes_list = ["continuous", "nominal", "ordinal"]
        invalid_dtypes_list = []
        for dt in parameters_dict["adjust_cols_dtypes"]:
            if dt not in valid_dtypes_list:
                invalid_dtypes_list.append(dt)
        
        if len(invalid_dtypes_list) > 0:
            log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###ERROR###: The following datatypes are not valid: {invalid_dtypes_list}.\n"
            helper.logger(log_file, log_string)
            sys.exit(1)


        # Check for missing values in columns for adjustment
        for var in parameters_dict["adjust_cols"]:
            if metadata_df[var].isna().any():
                n_na = metadata_df[metadata_df["MS.Batch"].isna()].shape[0]
                log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###ERROR###: Column {var} contains n={n_na} missing values.\n"
            helper.logger(log_file, log_string)
            sys.exit(1)



        # Check if datatype can be casted to float or string
        for var, dt in zip (
            parameters_dict["adjust_cols"], 
            parameters_dict["adjust_cols_dtypes"]):
            if dt == "continuous":
                try:
                    metadata_df[var] = metadata_df[var].astype(float)
                except:
                    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###ERROR###: Column {var} could not be casted to {dt}.\n"

            elif dt == "ordinal":
                try:
                    metadata_df[var] = metadata_df[var].astype(int)
                except:
                    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###ERROR###: Column {var} could not be casted to {dt}.\n"
            
            elif dt == "nominal":
                #nominal variable should be casted to string -> object in pandas
                try:
                    metadata_df[var] = metadata_df[var].astype(str)
                except:
                    log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"###ERROR###: Column {var} could not be casted to {dt}.\n"
            else:
                log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                +f"###ERROR###: Datatype not valid.\n"                

        # Check if requested model is implemented
        if parameters_dict["adjust_model"].lower() not in [
                "ols","ridge","lasso","elasticnet","huber","bayesian","loess"]:
                
                log_string += time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                +f"###ERROR###: Model {parameters_dict['adjust_model']} not implemented, no correction applied.\n"                

        if len(parameters_dict["adjust_fitsamples_meta"]) == 2:
            
            fit_col, fit_level = parameters_dict["adjust_fitsamples_meta"]
            fitsamples_idx_list = list(metadata_df[metadata_df[fit_col]==fit_level].index)
            info_fitsamples_str = \
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"Adjustment models will be fitted on {len(fitsamples_idx_list)} samples.\n"\
                    +f"Level {fit_level} in meta column {fit_col} will be used to adjust samples."  + "\n"
            log_string += info_fitsamples_str

            # Check whether the same levels are defined in samples

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
        
        if "ERROR" in log_string:
            # Exit modul because of error
            helper.logger(log_file,log_string)
            sys.exit(1)
        
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
                alphas_models = np.logspace(*parameters_dict["adjust_alphas_models"])
            )

    elif parameters_dict["correct_method"] == "adjust_confounders_batchwise":
        
        if parameters_dict["adjust_cols"] is None:
            warning_batch_level_str = \
                time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                +f"###ERROR###: No correction vriable defined \n"
            log_string += warning_batch_level_str
            sys.exit(1)

        fitsamples_idx_list = None
        if len(parameters_dict["adjust_fitsamples_meta"]) == 2:
            fit_col, fit_level = parameters_dict["adjust_fitsamples_meta"]
            fitsamples_idx_list = list(metadata_df[metadata_df[fit_col]==fit_level].index)
            info_fitsamples_str = \
                    time.strftime("%Y-%m-%d%t%H:%M:%S: ") \
                    +f"Correction model will be fitted on {len(fitsamples_idx_list)} samples.\n"\
                    +f"Level {fit_level} in meta column {fit_col} will be used to fit COMBAT Model."  + "\n"
            log_string += info_fitsamples_str

        data_processed_df, log_string, corr_dfs = adjust_confounders_batchwise(
            data_df,
            metadata_df, 
            model_string=parameters_dict["adjust_model"], 
            batch_level=parameters_dict["batch_cols_meta"], 
            random_state=seed, 
            covariates=parameters_dict["adjust_cols"], 
            covariates_dtypes_list=parameters_dict["adjust_cols_dtypes"],
            alpha_significance=0.05, alphas_models=[-3,3,7],
            fitsamples_idx_list=fitsamples_idx_list, correct_one_by_one=parameters_dict["adjust_batchwise_onebyone"])
    else:
        data_processed_df = data_df.copy()

    data_processed_df = data_processed_df.T

    helper.close_module(data_processed_df=data_processed_df, 
                        masked_df=masked_df, 
                        mask_description_df=mask_description_df, 
                        log_file=log_file, module_output=log_string,
                        out_dir_fp=out_dir_fp, outname=outname, mask=True)

# %%
