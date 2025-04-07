"""
Input: wide-format data
Output: Corrected wide-format data


- batchcorrection_combat: applies combat batch correction
- batchcorrection_harmony: applies harmony batch correction
- correct_for_confounders: drift correction per id for a given set of confounders which are statistically significant correlated to the id
- batchwise_drift_correction: applies drift correction per batch for a set of covariates



"""


#general
import pandas as pd
import numpy as np
import time
from tqdm import tqdm

#own scripts
from Python.Snippets import global_variables as global_variables # kept as place holder for logging, will be removed during programming of main 
from modules_python import module_helper

#visualization
import matplotlib.pyplot as plt

#statistics and modelling
import pycombat #CoAxLab
#from harmony import harmonize
from sklearn.preprocessing import StandardScaler
from scipy import stats
from sklearn import linear_model
from skmisc import loess
from sklearn import metrics





def batchcorrection_combat(profiles,correction_level,fitsamples=None,covariates=None,binlevel=False,bins=10,sample_column="sample"):
    """
    Apply combat batch correction, see https://github.com/CoAxLab/pycombat
    #TODO
    check whether covariates provided overlap with profiles.columns, if not inform the user
    
    author: BZ
    
    Parameters
    ----------
    profiles : pd.DataFrame
        features (e.g. precursors) as rows, samples as columns
    correction_level : str
        column level to correct by (e.g. "batch")
    fitsamples : list, optional
        samples to use for fitting the combat class (e.g. QC-samples only), must be column names of profiles. The default is None, fitting with all samples in profiles.
    covariates : list
        list of column levels with covariates to be corrected for
    binlevel : bool
        if True, bin the levels provided by correction_level
    bins : int
        number of bins for the binning of correction_level
    sample_column : str
        name of the column level with the sample names. The default is "sample".
        
    Returns
    -------
    pd.DataFrame, same shape as input profiles.

    """
    print("Batch correction using combat: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tBatch correction using combat\n")
        if fitsamples is not None:
            if type(fitsamples[0])!=str:
                fitsamplenames=fitsamples.get_level_values(sample_column)
            else:
                fitsamplenames=fitsamples
            logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t\t\t\tsamples used to fit batch correction:\n")
            logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t\t\t\t"+"; ".join(fitsamplenames)+"\n")
        logfile.close()
    #binning of level of interest
    if binlevel:
        counts,bins = np.histogram(profiles.columns.get_level_values(correction_level),bins=bins)
        bins[0]=bins[0]-0.01
        groups=pd.cut(profiles.columns.get_level_values(correction_level),bins,labels=range(1,len(bins)))
        profiles=module_helper.insert_column_level(df=profiles,new_level=groups.tolist(),level_name=correction_level)
    if fitsamples is None:
        fitsamples=profiles.columns
    combat=pycombat.Combat()
    #define covariates
    if covariates is not None:
        C_fit=profiles[fitsamples].columns.droplevel(list(set(profiles.columns.names) - set(covariates))).to_frame()
        combat.fit(Y=profiles.T.loc[fitsamples].values,b=profiles[fitsamples].columns.get_level_values(correction_level).tolist(),C=C_fit)
        C=profiles.columns.droplevel(list(set(profiles.columns.names) - set(covariates))).to_frame()
        data_corrected=combat.transform(Y=profiles.T.values,b=profiles.columns.get_level_values(correction_level).tolist(),C=C)
    else:
        combat.fit(Y=profiles.T.loc[fitsamples].values,b=profiles[fitsamples].columns.get_level_values(correction_level).tolist())
        data_corrected=combat.transform(Y=profiles.T.values,b=profiles.columns.get_level_values(correction_level).tolist())
    profiles_corrected=pd.DataFrame(data=data_corrected,columns=profiles.index,index=profiles.columns).T
    return(profiles_corrected)


def batchcorrection_harmony(profiles,correction_level,fitsamples=None,binlevel=False,bins=10):
    """
    Apply harmony batch correction, see https://pypi.org/project/harmony-pytorch/
    
    author: BZ
    
    Parameters
    ----------
    profiles : pd.DataFrame
        features (e.g. precursors) as rows, samples as columns
    correction_level : str
        column level to correct by (e.g. "batch")
    fitsamples : list, optional
        samples to use for fitting the combat class (e.g. QC-samples only), must be column names of profiles. The default is None, fitting with all samples in profiles.
    binlevel : bool
        if True, bin the levels provided by correction_level
    bins : int
        number of bins for the binning of correction_level

    Returns
    -------
    pd.DataFrame, same shape as input profiles.

    """
    
    from harmony import harmonize
    print("Batch correction using combat: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tBatch correction using combat\n")
        if fitsamples is not None:
            if type(fitsamples[0])!=str:
                fitsamplenames=fitsamples.get_level_values(0)
            else:
                fitsamplenames=fitsamples
            logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t\t\t\tsamples used to fit batch correction:\n")
            logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t\t\t\t"+"; ".join(fitsamplenames)+"\n")
        logfile.close()
    #binning of level of interest
    if binlevel:
        counts,bins = np.histogram(profiles.columns.get_level_values(correction_level),bins=bins)
        bins[0]=bins[0]-0.01
        groups=pd.cut(profiles.columns.get_level_values(correction_level),bins,labels=range(1,len(bins)))
        profiles=module_helper.insert_column_level(df=profiles,new_level=groups.tolist(),level_name=correction_level)
    if fitsamples is None:
        fitsamples=profiles.columns
    metadata_df=profiles.columns.to_frame()
    profiles_corrected = harmonize(profiles.T.values,metadata_df,batch_key=correction_level)
    profiles_corrected = pd.DataFrame(data=profiles_corrected,columns=profiles.index,index=profiles.columns).T
    return(profiles_corrected)



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


#
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
            r2_dfs.append(r2_df)
            #r2_df.to_csv(qc_folder + savefilename_precursors + appendix + "reprep_normalized_imputed,imputed_loesscorrected_on_all_for_SPBatch_"+str(int(spbatch))+".tsv",sep="\t")
            corrected_logprofiles.append(corr_logprofiles)
            #concatenate all corrected profiles
            logprofiles_corr=pd.concat(corrected_logprofiles,axis=1)
    r2_df_final=pd.concat(r2_dfs,axis=1)
    return([logprofiles_corr, r2_df_final])
    
    
