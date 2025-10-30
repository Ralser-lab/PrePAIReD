"""
Input: wide-format data
Output: imputed wide-format data

Missingness estimation functions
    - assign_missingness: determines whether no precursor per protein is measured
    - detection_probability_curve: Fits logistic regression to precursor detection probability and calculates a decision boundary for missingness classification
    - intensity_dependent_missingness_estimation: models Completeness over precursor intensity to classify precursor missingness into missing at random or not at random.
Imputation functions
    - graphbased_imputation: Graph-based approach to impute missing precursors.
    - knn_imputation: k-nearest neighbours imputation on peptide or protein matrix.
    - minval_imputation: minimum value imputation on peptide or protein matrix.
    - meanval_imputation: mean value imputation on peptide or protein matrix.
    - medianval_imputation: median value imputation on peptide or protein matrix.
    - bpca_fill: Bayesian PCA imputation on peptide or protein matrix.
Additional helper and logging functions
    - getIndexes: returns a list of positions where element exists in the dataframe
    - check_imputation_table_on_completeness: controls imputation_table on completeness for defined levels and values.
    - obtain_missingness_for_log: obtains list of tuples with positions of missing values and numbers of total and missing values in the given data_df.
    - update_mask: given a list of position tuples (locs) updates mask_df by adding mask_value.

"""


#general
import pandas as pd
import numpy as np
import re
import sys
import time
import json
import argparse
import os
import itertools
from collections import Counter

#statistics and modelling
import igraph
from sklearn.impute import KNNImputer
import statsmodels.api as sm
from scipy.optimize import curve_fit
from sklearn.preprocessing import SplineTransformer

#own scripts
#from modules_python import module_helper as helper
import module_helper as helper


########################################################################
# Missingness estimation functions
########################################################################

def assign_missingness(data_df,id_mapping_df,missingness_type="ProtAbsence",mask_df=None,precursor_col="Precursor.Id",protein_col="Protein.Group",on_filter_mask=True):
    """
    Determines whether no precursor of the proteins is measured and if so sets its value to zero.
    author: BZ
    

    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).
    id_mapping_df : pd.DataFrame
        dataframe mapping protein origin to precursors.
    missingness_type : str
        name of the missingness_type in missingness_types_df to use in returned missingness mask. The default is "ProtAbsence".
    mask_df : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None, generating the mask with the DIA-NN output file.
    precursor_col : string, optional
        name of the column in id_mapping_df containing precursor ids mapping to data_df and mask_df. The default is "Precursor.Id".
    protein_col : string, optional
        name of the column in id_mapping_df containing protein ids. The Default is "Protein.Group".
    on_filter_mask : bool, optional
        flag specifying whether to treat previously filtered precursors as present. The default is True. 

    Returns
    -------
    mask_df : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values corresponding to missingness as specified in missingness_types_df.

    """
    
    #remember previous sorting
    #orig_data_indices = data_df.index
    #orig_mask_indices = mask_df.index
    #copy original data for index operations
    mask_mi_df = mask_df.copy()
    data_mi_df = data_df.copy()
    data_mi_df.replace(np.nan,0,inplace=True)
    #make sure all measured values are above zero
    data_mi_df = data_mi_df.abs()
    #create multiindex for grouping
    data_mi_df.index = pd.MultiIndex.from_frame(id_mapping_df.set_index(precursor_col,drop=False).loc[data_mi_df.index,[precursor_col,protein_col]])
    mask_mi_df.index = pd.MultiIndex.from_frame(id_mapping_df.set_index(precursor_col,drop=False).loc[mask_mi_df.index,[precursor_col,protein_col]])
    #sort profiles by protein_col
    data_mi_df.sort_index(axis=0,level=protein_col,inplace=True)
    #create mask defining missing (False) or present (True) on protein-level
    if on_filter_mask: #filtered precursors count as presence of a protein
        missingness_mask_df = (mask_mi_df>=1).groupby(level=protein_col).sum()>0
        missingness_mask_df = missingness_mask_df.loc[data_mi_df.index.get_level_values(protein_col).unique(),data_mi_df.columns]
    else: #only consider precursors present after filtering
        missingness_mask_df = data_mi_df.groupby(level=protein_col).sum()>0
    
    #replace all missing values with -1
    missingness_df = data_mi_df.replace(0,-1)
    #multiply values with mask, resulting in 0 for values without any precursor from a given protein, -1 for values missing (likely) at random and >0 for measured values
    missingness_df = missingness_df.mul(missingness_mask_df,axis=0,level=protein_col)
    
    #revert labels such that False equals measured or missing and True equals missing without any precursor of given protein
    missingness_df = missingness_df==0
    
    #replace True with missingness label and False with zero
    missingness_df.replace(False,"",inplace=True)
    missingness_df.replace(True,missingness_type,inplace=True)
    return(missingness_df)
    

def detection_probability_curve(data_df, axis=1, boundary=0.5, detection_proportion_cap=1, method="logit-linear", n_knots=3):
    """
    Fits logistic regression to precursor detection probability and calculates a decision 
    boundary for missingness classification. Adapted from Li and Smyth, https://academic.oup.com/bioinformatics/article/39/5/btad200/7126416

    author: SB, BZ

    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).    
    axis : int
        1 for data_df precursors as rows, 0 for data_df with precursors as columns. The default is 1.
    boundary : float, optional
        Decision probability boundary for classification (default 0.5).
    detection_probability_cap : float
        A value 0 < detection_probability_cap <= 1 imposing a cap on the detection proportion for all precursors
    method : str
        A string specifying the method to fit the average precursor abundance to completeness. The deault is logit-linear.
    n_knots : int
        Number of knots used for data splines if method == logit-linear-splines. The default is 3.

    Returns
    -------
    decision_boundary : float
        Calculated abundance threshold for MNAR/MAR classification.
    result : statsmodels.GLMResults
        Logistic regression model fit results.
    """
    
    def func1_logitlinear(x,b0,b1):
        return b0 + b1 * x
    
    def func2_logitlinear_df(X,b):
        return np.sum(X*b,axis=1)
    
    def func3_capped_logitlinear(x,b0,b1):
        return detection_proportion_cap*np.exp(b0+b1*x)/(1+np.exp(b0+b1*x))

    def func4_logitlinear_biascorrected(X,b0,b1):
        return b0 + b1*X['Avg Log2 Abundance'] - (0.5*(b1**2)*X['Var Log2 Abundance'])
    
    #compute means, variances and missingness
    means = data_df.mean(axis=axis) 
    variances = data_df.var(axis=axis,ddof=1) 
    missingness = data_df.isnull().mean(axis=axis)
    
    #dataframe containing results
    res_df = pd.DataFrame({
        'Avg Log2 Abundance': means,
        'Var Log2 Abundance': variances,
        'Completeness': 1 - missingness,
        'Logit Completeness': np.log((1 - missingness)/missingness)
    }).replace([-np.inf,np.inf],np.nan).dropna(subset=["Avg Log2 Abundance","Completeness"])
    #sort res_df by abundance
    res_df.sort_values(by="Avg Log2 Abundance",inplace=True)
    
    if method=="logit-linear":
        #Linear model with binomial distribution for logit completeness
        X = sm.add_constant(res_df['Avg Log2 Abundance'])
        model = sm.GLM(res_df['Completeness'], X, family=sm.families.Binomial())
        result = model.fit()
        intercept, coef = result.params
        res_df["Completeness pred"] = 1/(1+np.exp(-func1_logitlinear(res_df["Avg Log2 Abundance"],intercept,coef)))
        #decision_boundary = (np.log(boundary / (1 - boundary)) - intercept) / coef 
    elif method=="logit-linear-fit":
        popt, pcov = curve_fit(func1_logitlinear,res_df.dropna(subset='Logit Completeness')['Avg Log2 Abundance'],res_df.dropna(subset='Logit Completeness')['Logit Completeness'],method="trf",loss="cauchy")
        res_df["Completeness pred"] = 1/(1+np.exp(-func1_logitlinear(res_df['Avg Log2 Abundance'],*popt)))
        #decision_boundary =  (np.log(boundary / (1 - boundary)) - popt[0]) / popt[1]
    elif method=="logit-linear-splines":
        indices = res_df.dropna(subset=["Avg Log2 Abundance","Completeness","Logit Completeness"]).index
        X = sm.add_constant(res_df.loc[indices,'Avg Log2 Abundance'])
        #model = make_pipeline(SplineTransformer(n_knots=n_knots, degree=1), sm.GLM(family=sm.families.Binomial()))
        spline = SplineTransformer(n_knots=n_knots, degree=1)
        splines = spline.fit_transform(res_df.loc[indices,'Avg Log2 Abundance'].to_numpy().reshape(-1,1),res_df.loc[indices,"Logit Completeness"])
        X = sm.add_constant(splines)
        model = sm.GLM(res_df.loc[indices,'Completeness'], X, family=sm.families.Binomial())
        result = model.fit()
        popt = result.params
        res_df.loc[indices,"Completeness pred"] = 1/(1+np.exp(-func2_logitlinear_df(X=X,b=popt.to_numpy().reshape(1,-1))))
    elif method=="logit-linear-capped":
        popt, pcov = curve_fit(func3_capped_logitlinear,res_df['Avg Log2 Abundance'],res_df['Completeness'],bounds=([-np.inf,0],[np.inf,1]))
        res_df["Completeness pred"] = func3_capped_logitlinear(res_df['Avg Log2 Abundance'],*popt)
    elif method=="logit-linear-biascorrected":
        indices = res_df.dropna().index
        popt, pcov = curve_fit(func4_logitlinear_biascorrected,res_df.loc[indices],res_df.loc[indices,'Logit Completeness'])
        #res_df.loc[indices,"Completeness pred"] = 1/(1+np.exp(-func4_logitlinear_biascorrected(res_df.loc[indices],*popt)))
        #fit too fuzzy, therefore logit-linear fit on predicted values
        logit_completeness = func4_logitlinear_biascorrected(res_df.dropna(),*popt)
        popt, pcov = curve_fit(func1_logitlinear,res_df.loc[indices,"Avg Log2 Abundance"],logit_completeness)
        res_df.loc[indices,"Completeness pred"] = 1/(1+np.exp(-func1_logitlinear(res_df.loc[indices,"Avg Log2 Abundance"],*popt)))
    elif method=="response-binned-splines": #from SB
        # Spline the response variable (Detection Probability) using pd.qcut()
        res_df['response_bin'] = pd.qcut(res_df['Completeness'], q=100, duplicates='drop')
        # Compute mean abundance and mean detection probability per bin
        grouped = res_df.groupby('response_bin').agg(abundance_mean=('Avg Log2 Abundance', 'median'),detection_mean=('Completeness', 'median')).dropna()
        # Prepare data for GLM
        X = sm.add_constant(grouped['abundance_mean'])
        y = grouped['detection_mean']
        # Fit binomial GLM (using proportions directly with freq_weights)
        model = sm.GLM(y, X, family=sm.families.Binomial())
        result = model.fit()
        # Decision boundary calculation
        intercept, coef = result.params
        res_df["Completeness pred"] = 1/(1+np.exp(-func1_logitlinear(res_df["Avg Log2 Abundance"],intercept,coef)))
    elif method=="logit-linear-capped-two-sides":
        #include asymptotic cap on lower completeness and/or weighting by thresholding lower completeness values
        pass #TODO 
        
    decision_boundary = res_df.loc[(res_df["Completeness pred"].dropna()>=boundary).idxmax(),"Avg Log2 Abundance"]
    
    return round(decision_boundary, 2), res_df


def intensity_dependent_missingness_estimation(data_df, axis=1, boundary=0.5, detection_proportion_cap=1, method="logit-linear", n_knots=3):
    """
    Integrated pipeline for preprocessing data, classifying precursors, and performing 
    mixed imputation (minimum and KNN).

    author: SB, BZ

    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities (e.g. Precursor.Quantity from DIA-NN output).    
    axis : int
        1 for data_df precursors as rows, 0 for data_df with precursors as columns. The default is 1.
    boundary : float, optional
        Decision probability boundary for classification, i. e. Completeness threshold to obtain abundance boundary. The default is 0.5.
    detection_probability_cap : float
        A value 0 < detection_probability_cap <= 1 imposing a cap on the detection proportion for all precursors
    method : str
        A string specifying the method to fit the average precursor abundance to completeness. The deault is logit-linear.
    n_knots : int
        Number of knots used for data splines if method == logit-linear-splines. The default is 3.


    Returns
    -------
    mnar_mask df : pd.DataFrame
        Dataframe specifying missingness for every precursor in every sample.
    regression_df : pd.DataFrame
        Dataframe with average precursor abundance, variance, completeness and predicted completeness by fit function.
        
    """
    
    # Fit Logistic Regression
    decision_boundary, regression_df = detection_probability_curve(data_df, axis=axis, boundary=boundary, detection_proportion_cap=detection_proportion_cap, method=method, n_knots=n_knots)
    
    mar_int = 2
    mnar_int = 3
    
    #initialize mask distinguishing MAR from MNAR, preserving measured values as zeros
    missingness_mask_df = data_df.isnull()
    missingness_mask_df.replace({False: 0, True: 1},inplace=True)
    #classify missingness type for every precursor as MAR if below decision_boundary and MNAR if above decision_boundary
    mars = (data_df.mean(axis=axis) < decision_boundary)
    mars = mars[mars].index
    mnars = (data_df.mean(axis=axis) >= decision_boundary)
    mnars = mnars[mnars].index
    missingness_mask_df.loc[mars] = missingness_mask_df.loc[mars] * mar_int
    missingness_mask_df.loc[mnars] = missingness_mask_df.loc[mnars] * mnar_int
    #reset values which are missing but missingness could not be classified back to zero
    missingness_mask_df.replace(1,0,inplace=True)
    #missingness_mask_df.replace(0,np.nan,inplace=True)
    #numbers to strings
    missingness_mask_df.replace({np.nan:"", 0:"", mar_int:"MAR", mnar_int:"MNAR"},inplace=True)
    #missingness_mask_df.replace([0,np.nan],"",inplace=True)
    #missingness_mask_df.replace(mar_int,"MAR",inplace=True)
    #missingness_mask_df.replace(mnar_int,"MNAR",inplace=True)
    
    return missingness_mask_df, regression_df, decision_boundary


########################################################################
# Imputation functions
########################################################################

def graphbased_imputation(data_df,id_mapping_df,protein_col="Protein.Group",precursor_col="Precursor.Id",impute_modifications=True,estimate_threshold=1,weighting=False):
    """
    Graph-based approach to impute missing precursors. This approach assumes that the ratio of any two precursors from a given protein solely relies on physico-chemical properties and digestion and therefore is equal in all samples/runs.
    Creates for every protein group a directed, weighted graph with precursors as nodes and median ratio of any two precursors (nodes) over all samples/fraction as edge weights between nodes.
    These ratios are then multiplied with present precursors over samples, median over them is used as intensity of missing precursor. 
    This approach may not impute all missing precursors, as advantage it uses global information over all samples about precursor ratios, but imputes missing values only from precursor intensities within a sample times global ratios.
    Imputed values are marked in mask_df as a second number after decimal point in the matrix, preserving the information if the precursor was previously filtered or not measured at all
    
    author: BZ
    
    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    id_mapping_df : pd.DataFrame
        dataframe mapping protein origin to precursors.
    protein_col : string, optional
        name of the column in id_mapping_df containing protein ids. The Default is "Protein.Group".
    precursor_col : string, optional
        name of the column in id_mapping_df containing precursor ids mapping to data_df and mask_df. The default is "Precursor.Id".
    impute_modifications : bool, optional
        Switch to True to also impute precursors with modifications. In fractionation experiments the ratio between modified precursors is not necessarily equal over fractions(=samples). The default is True.
    estimate_threshold : int, optional
        Minimum number of annotated edges (known ratios between two precursors) to impute a missing precursor value. The default is 1.
    weighting : bool, optional
        If True calculates the weighted median between estimated imputed values for a precursor weighted on how often the ratio between two precursors was used to calculate edge weights. The default is False.
        
    Returns
    -------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities with imputed values.

    """
   
    #for every protein group
    for pg in id_mapping_df[protein_col].unique():
        #create graph for protein with precursors as nodes and ratios as edges
        temp_df = data_df[id_mapping_df.loc[protein_col]==pg,precursor_col]
        #exclude columns without any precursor intensities for given protein group
        temp_df.dropna(how="all",axis=1,inplace=True)
        #exclude precursors with modification sites
        if impute_modifications==False:
            modi_precs={x if "UniMod:" in x else np.nan for x in temp_df.index}
            if np.nan in modi_precs:
                modi_precs.remove(np.nan)
            if modi_precs!=set():
                #include non-modified sequence
                modi_precs=modi_precs | {re.sub("\(.*\)","",prec) for prec in modi_precs}
                temp_df.drop(temp_df[temp_df.index.isin(modi_precs)].index,inplace=True)
        if temp_df.empty:
            continue
        #edges=permutations(range(len(temp_df.index)),2)
        #edge_weights=list(map(lambda x: np.nanmedian(temp_df.loc[x[0]]/temp_df.loc[x[1]]),permutations(temp_df.index,2)))
        #ratio_weights tell how often the ratio between two given precursors was used to calculate edge weight. Edges with higher ratio_weight have higher impact in calculation of median for the missing value
        g=igraph.Graph(directed=True)
        g.add_vertices(range(len(temp_df.index)))
        g.vs["name"]=temp_df.index
        g.add_edges(itertools.permutations(range(len(temp_df.index)),2))
        #g.es["weight"]=list(map(lambda x: np.nanmedian(temp_df.loc[x[1]].astype(float)/temp_df.loc[x[0]].astype(float)),permutations(temp_df.index,2)))
        g.es["weight"]=list(map(lambda x: np.nanmedian(temp_df.loc[x[1]]/temp_df.loc[x[0]]),itertools.permutations(temp_df.index,2)))
        if weighting:
            g.es["ratio_weight"]=list(map(lambda x: len((temp_df.loc[x[1]]/temp_df.loc[x[0]]).dropna()),itertools.permutations(temp_df.index,2)))
        g.es["name"]=list(itertools.permutations(temp_df.index,2))
        '''
        #plot
        g.es["label"]=g.es["weight"]
        fig, ax = plt.subplots()
        plt.axis("off")
        igraph.plot(g,target=ax,bbox=(300,300),vertex_color=["#4E79A7"]*g.vcount(),vertex_size=15,autocurve=False,vertex_label=[x[0] for x in g.vs["name"]])
        '''
        #if no edges in graph, skip
        #...
        #impute to every sample/fraction missing precursors with graph edge weights
        for col in temp_df.columns[temp_df.isna().any()]: #sample/fraction to impute precursors in
            #df[nan-precursor,sample]=mean(rel_abundance)*present_precursors #mean(rel_abundance) is list of lists, for every precursor over all samples
            for precursor in temp_df[temp_df[col].isna()].index:
                #estimates of missing value from its neighbour (predecessor) nodes
                estimates=[g.es.find(name=(x["name"],precursor))["weight"]*temp_df.loc[x["name"],col] for x in g.vs.find(name=precursor).predecessors()]
                if weighting:
                    ratio_weights=[g.es.find(name=(x["name"],precursor))["ratio_weight"]*temp_df.loc[x["name"],col] for x in g.vs.find(name=precursor).predecessors()]
                #calculate weighted median to estimate missing value, if at least estimate_threshold estimates for missing value are given
                if np.sum(~np.isnan(estimates))>=estimate_threshold:
                    if weighting:
                        data_df.loc[precursor,col]=helper.weighted_nanmedian(estimates,ratio_weights)
                    else:
                        data_df.loc[precursor,col]=np.nanmedian(estimates)
                #[g.es.find(name=(x["name"],precursor))["weight"]*temp_df.loc[x["name"],col] for x in g.vs.find(name=precursor).predecessors()]
                #[g.es.find(name=(x["name"],precursor))["weight"] for x in g.vs.find(name=precursor).predecessors()]
                #[temp_df.loc[x["name"],col] for x in g.vs.find(name=precursor).predecessors()]
    return(data_df)



def knn_imputation(data_df,k=5,weights="distance",metric="nan_euclidean",transposed=True):
    """
    k-nearest neighbours imputation on peptide or protein profiles.

    author: BZ
    
    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    k : int, optional
        number of neighbours to take into account for kNN Imputation. The default is 5.
    weights : str, optional
        Weight function used in prediction as passed to KNNImputer. The default is distance.
    metric : str, optional
        Distance metric for searching neighbors as passed to KNNImputer. The default is "nan_euclidean".
    transposed : bool, optional
        If True, tranposes the data matrix prior to imputation and impute as described by Troyanskaja et al.https://doi.org/10.1093/bioinformatics/17.6.520 , else imputes missing values of a protein based on values of k neighbouring samples

    Returns
    -------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.


    """
    #knn imputation
    imputer = KNNImputer(n_neighbors=k,weights=weights,keep_empty_features=True, metric=metric)
    indices=data_df.index
    columns=data_df.columns
    if transposed:
        imputed_values = imputer.fit_transform(data_df)
        data_df = pd.DataFrame(index=indices,columns=columns,data=imputed_values)
    else:
        imputed_values = imputer.fit_transform(data_df.T)
        data_df = pd.DataFrame(index=columns,columns=indices,data=imputed_values).T
    return(data_df)
   

def minval_imputation(data_df,axis=1,add_noise=True):
    """
    minimum value imputation on peptide or protein profiles. Each missing value is imputed by the minimum of the row, if add_noise=True divided by the square root of 2 times a random number between 0.75 and 1.25.
    
    author: BZ
    
    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    axis : int
        1 for data_df precursors as rows, 0 for data_df with precursors as columns. The default is 1.
    add_noise : bool, optional
        divides the minimum value by the square root of 2 (1.41) and multiplies it by a random uniform value between 0.75 and 1.25, thus selecting a value slightly below the minimum. The default is True.
    Returns
    -------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    """
    
    #impute by minimum value per ID divided by sqrt(2) times random number between 0.75 and 1.25
    if add_noise:
        if axis==1:
            for row in data_df.index:
                data_df.loc[row]=data_df.loc[row].apply(lambda x: data_df.loc[row].min()/np.sqrt(2)*np.random.uniform(0.75,1.25) if pd.isna(x) else x)
        elif axis==0:
            for col in data_df.columns:
                data_df[col]=data_df[col].apply(lambda x: data_df[col].min()/np.sqrt(2)*np.random.uniform(0.75,1.25) if pd.isna(x) else x)
    else:
        #simple minimum value imputation
        if axis==1:
            for row in data_df.index:
                data_df.loc[row] = data_df.loc[row].fillna(data_df.loc[row].min())
                #data_df.fillna(data_df.min(axis=axis),axis=axis,inplace=True) #NotImplementedError for axis=1
        else:
            data_df.fillna(data_df.min(axis=axis),axis=axis,inplace=True)
    return(data_df)
    

def meanval_imputation(data_df,axis=1,add_noise=True):
    """
    mean value imputation on peptide or protein profiles. Each missing value is imputed by the mean of the row, if add_noise=True with a value between mean - standard deviation and mean + standard deviation.
    
    author: BZ
    
    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    axis : int
        1 for data_df precursors as rows, 0 for data_df with precursors as columns. The default is 1.
    add_noise : bool, optional
        adds to the mean value a random normal value between -standard deviation and +standard deviation, thus selecting a value around the mean. The default is True.

    Returns
    -------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.

    """
    #impute by mean value per ID
    if add_noise:
        if axis==1:
            for row in data_df.index:
                data_df.loc[row]=data_df.loc[row].apply(lambda x: data_df.loc[row].mean() + np.random.normal(-data_df.loc[row].std(),data_df.loc[row].std()) if pd.isna(x) else x)
        elif axis==0:
            for col in data_df.columns:
                data_df[col]=data_df[col].apply(lambda x: data_df[col].mean() + np.random.normal(-data_df[col].std(),data_df[col].std()) if pd.isna(x) else x)
    else:
        #simple mean value imputation
        if axis==1:
            for row in data_df.index:
                data_df.loc[row] = data_df.loc[row].fillna(data_df.loc[row].mean())
                #data_df.fillna(data_df.mean(axis=axis),axis=axis,inplace=True) #NotImplementedError for axis=1
        else:
            data_df.fillna(data_df.mean(axis=axis),axis=axis,inplace=True)
    return(data_df)


def medianval_imputation(data_df,axis=1,add_noise=True):
    """
    median value imputation on peptide or protein profiles. Each missing value is imputed by the median of the row, if add_noise=True with a value between median - MAD and median + MAD.
    
    author: BZ
    
    Parameters
    ----------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    axis : int
        1 for data_df precursors as rows, 0 for data_df with precursors as columns. The default is 1.
    add_noise : bool, optional
        adds to the median value a random normal value between -MAD and +MAD, thus selecting a value around the median. The default is True.

    Returns
    -------
    data_df : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.

    """
    #impute by median value per ID
    if add_noise:
        if axis==1:
            for row in data_df.index:
                data_df.loc[row]=data_df.loc[row].apply(lambda x: data_df.loc[row].median() + np.random.normal(-helper.mad(data_df.loc[row]),helper.mad(data_df.loc[row])) if pd.isna(x) else x)
        elif axis==0:
            for col in data_df.columns:
                data_df[col]=data_df[col].apply(lambda x: data_df[col].median() + np.random.normal(-helper.mad(data_df[col]),helper.mad(data_df[col])) if pd.isna(x) else x)
    else:
        #simple median value imputation
        if axis==1:
            for row in data_df.index:
                data_df.loc[row] = data_df.loc[row].fillna(data_df.loc[row].median())
                #data_df.fillna(data_df.median(axis=axis),axis=axis,inplace=True) #NotImplementedError for axis=1
        else:
            data_df.fillna(data_df.median(axis=axis),axis=axis,inplace=True)
    return(data_df)


def bpca_fill(dX, max_iter = 200):
     """
    Bayesian PCA imputation on peptide or protein profiles as described by Oba et al. https://doi.org/10.1093/bioinformatics/btg287
    downloaded from https://github.com/unifey/bPCAfill, 18.03.2024 17:08, credit to unifey

    Parameters
    ----------
    dX : pd.DataFrame
        rows equal samples, columns ids (precursors or proteins).
    max_iter : int, optional
        maximum number of iterations. The default is 200.

    Returns
    -------
    dX_filled : pd.DataFrame
        rows equal samples, columns ids (precursors or proteins).
    """
    y = dX.copy() # Just following the conventions used in the original paper by Oba et al. (2003)
    cols = y.columns #.tolist()
    indices = y.index

    #--- Normalize `y` and convert to matrix format 
    means = y.mean(); sd = y.std(); y = (y-means)/sd
    y = np.matrix(y).astype(float) 

    #--- Initialization for Bayesian PCA imputation
    N,d = y.shape; q = d-1
    yest = np.copy(y); yest[np.isnan(yest)] = 0 

    missidx = {} # dictionary that contains NaN positions (row,col) of dX; {row : list of NaN col in row}
    bad = np.where(np.isnan(y))
    badrows = bad[0]; badcols = bad[1]
    for row in badrows: missidx[row] = []
    for link in range(len(badrows)): missidx[badrows[link]].append(badcols[link])

    nomissidx = {} # dictionary of indices complement to missidx
    good = np.where(~np.isnan(y))
    goodrows = good[0]; goodcols = good[1]
    for row in goodrows: nomissidx[row] = []
    for link in range(len(goodrows)): nomissidx[goodrows[link]].append(goodcols[link])

    ##--- Removing redundancies in row entries in `badrows` and `goodrows`
    gmiss = list(set(badrows)) # There are NaN(s) in these rows
    gnomiss = list(set(goodrows)) # None of these rows are completely NaN.

    covy = np.cov(yest.T) # covariance matrix of dX
    U, S, V = np.linalg.svd(np.matrix(covy)) # Singular Value Decomposition of the covariance matrix
    U = (U.T[0:q]).T;         S = S[0:q]*np.eye(q);           V = (V.T[0:q]).T

    mu = np.nanmean(np.copy(y), 0)
    W = U*np.sqrt(S)
    taumax = 1e20; taumin = 1e-20
    tau = 1/ (np.trace(covy)-np.trace(S));    tau = np.amax([np.amin([tau,taumax]),taumin])

    galpha0 = 1e-10;    balpha0 = 1;    alpha = (2*galpha0 + d)/(tau*np.diag(W.T*W)+2*galpha0/balpha0)
    gmu0  = 0.001;    btau0 = 1;    gtau0 = 1e-10;    SigW = np.eye(q)

    #--- End of bPCA initialization
    
    ### Outcome: N, d, q, yest, missidx, nomissidx, gmiss, gnomiss, 
    ###          mu, W, tau, galpha0, balpha0, alpha, gmu0, btau0, gtau0, SigW

    #--- Iteration

    ##--- Initialize `tau_old` and fill NaNs in `y` with `sys.maxsize`
    tau_old = 1000
    y = np.nan_to_num(y, nan = sys.maxsize)

    for epoch in range(max_iter):

        #--- Begin BPCA_dostep(): batch 1 step of Bayesian PCA EM algorithm
        Rx = np.eye(q)+tau*W.T*W+SigW;            Rxinv = np.linalg.inv(Rx)
        idx = gnomiss; n = len(idx)                  
        dy = y[idx,:] - np.tile(mu,(n,1));      x = tau * Rxinv * W.T * dy.T

        Td = dy.T*x.T;                            trS = np.sum(np.multiply(dy,dy))
        for n in range(len(gmiss)):
            i = gmiss[n]
            dyo = np.copy(y)[i,nomissidx[i]] - mu[nomissidx[i]]
            Wm = W[missidx[i],:];                                  Wo = W[nomissidx[i],:]
            Rxinv = np.linalg.inv( Rx - tau*Wm.T*Wm );             ex = tau * Wo.T * np.matrix(dyo).T;   x = Rxinv * ex
            dym = Wm * x;                                          dy = np.copy(y)[i,:]
            dy[nomissidx[i]] = dyo
            dy[missidx[i]] = dym.T
            yest[i,:] = dy + mu
            Td = Td + np.matrix(dy).T*x.T;                            Td[missidx[i],:] = Td[missidx[i],:] + Wm * Rxinv
            trS = trS + dy*np.matrix(dy).T +  len(missidx[i])/tau + np.trace( Wm * Rxinv * Wm.T )

        Td = Td/N;                trS = trS/N;

        ###--- Update outcome: (N, d, q,) yest(updated in n-loop above), (missidx, nomissidx, gmiss, gnomiss,)
        ###                 (mu,) W, tau, (galpha0, balpha0,) alpha, (gmu0, btau0, gtau0,) SigW

        Rxinv = np.linalg.inv(Rx); 
        Dw = Rxinv + tau*Td.T*W*Rxinv + np.diag(alpha)/N;    Dwinv = np.linalg.inv(Dw);

        W = Td * Dwinv;
        tau = (d+2*gtau0/N)/(trS-np.trace(Td.T*W)  + (mu*np.matrix(mu).T*gmu0+2*gtau0/btau0)/N)[0,0];
        SigW = Dwinv*(d/N);
        alpha = (2*galpha0 + d)/ (tau*np.diag(W.T*W)+np.diag(SigW)+2*galpha0/balpha0)

        dtau = np.log10(tau) - np.log10(tau_old)
        if np.abs(dtau) < 1e-4:  break;

        # Update `tau_old`
        tau_old = tau

    #--- Create imputed dataframe `dX_filled`
    
    dXb = pd.DataFrame(yest)
    dXb.columns = cols
    dXb.index = indices
    dXb = (dXb*sd)+means

    dXo = dX.copy()
    missing = dXo.isna()
    dXo[missing] = 0
    dXb[~missing] = 0

    dX_filled = dXo + dXb
    
    return dX_filled


########################################################################
# Additional helper and logging functions
########################################################################

def getIndexes(dfObj, value):
    """
    Returns a list of positions where element exists in the dataframe.
    copied from https://www.geeksforgeeks.org/find-location-of-an-element-in-pandas-dataframe-in-python/ , access 10.02.2023, 11:44
    
    Parameters
    ----------
    dfObj : pd.DataFrame
        dataframe in which positions for given value are searched.
    value : str, int or float
        element which positions are of interest.

    Returns
    -------
    listOfPos : list
        All positions where dfObj.loc[listOfPos]==value.

    """
    # Empty list
    listOfPos = []
    # isin() method will return a dataframe with
    # boolean values, True at the positions   
    # where element exists
    result = dfObj.isin([value])
    # any() method will return
    # a boolean series
    seriesObj = result.any()
    # Get list of column names where
    # element exists
    columnNames = list(seriesObj[seriesObj == True].index)
    # Iterate over the list of columns and
    # extract the row index where element exists
    for col in columnNames:
        rows = list(result[col][result[col] == True].index)
        for row in rows:
            listOfPos.append((row, col))
    # This list contains a list tuples with
    # the index of element in the dataframe
    return listOfPos


def check_imputation_table_on_completeness(impu_table_df,metadata_df,data_df=None,metadata_sample_column="Sample.Name"):
    if data_df is not None:
        samples = data_df.columns
    else:
        samples = metadata_df[metadata_sample_column]
    sample_hierarchies = impu_table_df.columns[:impu_table_df.columns.tolist().index("precursor_hierarchy_values")] 
    #generate list with dictionaries through hierarchies from metadata
    all_sample_hierarchy_values_dicts_list = list()
    sample_hierarchies_stack = np.copy(sample_hierarchies)
    while len(sample_hierarchies_stack)>1:
        sample_hierarchy = sample_hierarchies_stack[0]
        lower_sample_hierarchy = sample_hierarchies_stack[1]
        sample_hierarchies_stack = sample_hierarchies_stack[1:]
        sample_hierarchy_values = metadata_df[metadata_df[metadata_sample_column].isin(samples)][sample_hierarchy].unique()
        sample_hierarchy_values = np.sort(sample_hierarchy_values)
        sample_hierarchy_value_dict = dict()
        for sample_hierarchy_value in sample_hierarchy_values:
            lower_sample_hierarchy_values = metadata_df.loc[(metadata_df[metadata_sample_column].isin(samples)) & (metadata_df[sample_hierarchy]==sample_hierarchy_value)][lower_sample_hierarchy].unique()
            lower_sample_hierarchy_values = np.sort(lower_sample_hierarchy_values)
            sample_hierarchy_value_dict[sample_hierarchy_value] = lower_sample_hierarchy_values
        all_sample_hierarchy_values_dicts_list.append(sample_hierarchy_value_dict)
    all_sample_hierarchy_values_dicts_list = [{"":list(all_sample_hierarchy_values_dicts_list[0].keys())}] + all_sample_hierarchy_values_dicts_list
    #compare list with dictionaries through hierarchies with values in imputation table
    for h,sample_hierarchy in enumerate(sample_hierarchies):
        if h == 0:
            upper_hierarchy = sample_hierarchy
        #lower_hierarchy #TODO
        


def obtain_missingness_for_log(data_df):
    """
    Obtains list of tuples with positions of missing values and numbers of total and missing values in the given data_df.

    Parameters
    ----------
    data_df : pd.DataFrame
        Dataframe to identify missing values in.

    Returns
    -------
    total_num_values : int
        number of values in data_df
    num_missing_values : int
        number of missing values in data_df
    missing_val_locs : list
        list of tuples with positions of missing values in data_df

    """
    total_num_values = data_df.shape[0] * data_df.shape[1]
    num_missing_values = pd.isna(data_df).sum().sum()
    #remember which values were missing prior imputation
    missing_val_locs = getIndexes(data_df.isna(),True)
    return([total_num_values,num_missing_values,missing_val_locs])


def update_mask(mask_df, locs, mask_value):
    """
    Given a list of position tuples (locs) updates mask_df by adding mask_value.

    Parameters
    ----------
    mask_df : pd.DataFrame
        Dataframe with locs being a subset of its positions.
    locs : list
        list of tuples, every tuple is a index and a column in mask_df.
    mask_value : float, int or str
        A value that can be added to the values in mask_df, if mask_df contains only numeric values, mask_value must be numeric.

    Returns
    -------
    mask_df : pd.DataFrame

    """
    for precursor,column in locs:
        mask_df.loc[precursor,column] += mask_value
    return mask_df

    

if __name__ == "__main__":
    
    data_df, parameters_dict, id_mapping_df, seed, module_output, metadata_df, masked_df, mask_description_df, out_dir_fp, outname, log_file, counter = helper.init_module(module = "Imputation", 
                                                                                                                                                                           load_metadata=True, 
                                                                                                                                                                           load_id_mapping=True, 
                                                                                                                                                                           load_mask=True, 
                                                                                                                                                                           return_data_fp_only=False)
 
        
    #temporary for debugging while load_parameters is not ready, careful with eval if variables are already defined, quick and dirty
    '''
    meta_df_fp = "../Metadata.tsv"
    data_df_fp = "../01_Filtering/data_1.tsv"
    out_dir_fp = "../03_Imputation/"
    log_file = "../log.txt"
    parameters_dict_raw = pd.read_csv("../03_Imputation/parameters.txt",sep="\t",index_col=0,header=None)[1].to_dict()
    parameters_dict = dict()
    for x in parameters_dict_raw.items():
        try: #if type != str
            parameters_dict[x[0]] = eval(x[1])
        except (SyntaxError, NameError): #if type==str
            parameters_dict[x[0]] = x[1].split(" #")[0]
    
    # Load data
    data_df = pd.read_csv(data_df_fp,sep="\t",index_col=0)
    
    # Load ID mapping from Precursors to Proteins
    id_mapping_df = pd.read_csv(out_dir_fp + "../01_Filtering/id_mapping.tsv", sep="\t")
    
    # Load mask file
    mask_df = pd.read_csv(data_df_fp[:-4] + "_mask.tsv",sep="\t",index_col=0)
    mask_df = mask_df.astype(float)
    mask_description_df = pd.read_csv(data_df_fp[:-4] + "_mask_description.tsv",sep="\t",index_col=0)
    
    # Load metadata_df
    metadata_df = pd.read_csv(meta_df_fp,sep="\t")
    
    #metadata adaption to fit hierarchical imputation
    metadata_df["Batch"]=metadata_df["MS.Batch"]
    metadata_df["Batch_Plate"]=metadata_df["MS.Batch"].astype(str) + "_" + metadata_df["Plate"].astype(str)
    metadata_df["Batch_Plate_Row"] = metadata_df["MS.Batch"].astype(str) + "_" + metadata_df["Plate"].astype(str) + "_" + metadata_df["Plate.Row"].astype(str)
    '''
    
    # Define filepath of imputation procedures table #May be passed from wrapper? #TODO
    imputation_table_fp = out_dir_fp + "imputation_procedures.tsv"
    
    # Log parameters
    log = time.strftime("%Y-%m-%d%t%H:%M:%S")+"\t"+"Imputation\n"
    log += "\nParameters\n------------\n"
    log += "\t" + pd.DataFrame(parameters_dict.items(),columns=["Parameter","Value"]).to_string(index=False).replace("\n","\n\t") 
    log += "\n------------\n"  
    helper.logger(log_file, module_output=log)
    #log = ""
    
    
    # replace zeros with NA
    if parameters_dict["replace_zeros_with_nan"]:
        data_df.replace(0,np.nan,inplace=True)
    
    # only proceed if any value is missing
    if data_df.isnull().values.any():
        #create impu_table template
        if parameters_dict["hierarchical_imputation"]:
            if os.path.exists(imputation_table_fp) and not parameters_dict["overwrite_imputation_procedures"]:
                # load imputation table
                impu_table_df = pd.read_csv(imputation_table_fp,sep="\t")
                # check imputation table on completeness such that every lower hierarchy is defined for every upper hierarchy
                #check_imputation_table_on_completeness() #TODO
            else:
                # if none existent, create table with imputation methods
                sample_hierarchies = parameters_dict["imputation_sample_hierarchies"].split(";")
                #precursor_hierarchy_values = [float(x) for x in parameters_dict["missingness_imputation_thresholds"].split(";")]
                impu_table_df = helper.create_imputation_table(metadata_df=metadata_df,
                                                        sample_hierarchies=sample_hierarchies,
                                                        missingness_thresholds=parameters_dict["missingness_imputation_thresholds"],
                                                        data_df=data_df,
                                                        metadata_sample_column=parameters_dict["metadata_sample_column"],
                                                        default_impu_method_ProtAbsence=parameters_dict["imputation_method_ProtAbsence"],
                                                        default_impu_method_MAR=parameters_dict["imputation_method_MAR"],
                                                        default_impu_method_MNAR=parameters_dict["imputation_method_MNAR"],
                                                        default_impu_method_ProtAbsenceMAR=parameters_dict["imputation_method_ProtAbsenceMAR"],
                                                        default_impu_method_ProtAbsenceMNAR=parameters_dict["imputation_method_ProtAbsenceMNAR"],
                                                        default_knn_k=parameters_dict["knn_k"])
                impu_table_df.to_csv(imputation_table_fp,sep="\t",index=False)

            '''
            #extract missingness types 
            missingness_types_all = impu_table_df["missingness_type"].drop_duplicates().tolist()
            if parameters_dict["remember_missingness_type_in_mask"]:
                missingness_types_df = pd.DataFrame(columns=["filter","value"])
                missingness_types_df["filter"] = missingness_types_all
                missingness_types_df["value"] = [x/10 for x in range(1,(len(missingness_types_all)+1))]
                mask_description_df = pd.concat([mask_description_df,missingness_types_df.set_index("filter")],axis=0)
            '''
            
        # Process data
        try:
            # create a mask to store missingness type (thus imputation approach neccessary), 0 equals no imputation applied
            if parameters_dict["eval_missingness_on_protein_presence"]:
                #check that parameters_dict["missingness_type_protein_absence"] is in missingness types
                if parameters_dict["missingness_type_protein_absence"] not in missingness_types_all and parameters_dict["hierarchical_imputation"]:
                    log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\t###WARNING### Missingness evaluation on protein presence as defined by missingness_type_protein_absence in the parameters does not match any entry in imputation_table['missingness_type'].\n\t\t\t\tMissingness evaluation on protein presence skipped.\n"
                else:
                    # remember in mask missingness type (protein present vs. protein not present or MAR vs. MAR/MNAR)
                    protmissing_mask_df = assign_missingness(data_df=data_df,
                                                            id_mapping_df=id_mapping_df,
                                                            missingness_type=parameters_dict["missingness_type_protein_absence"],
                                                            mask_df=mask_df,
                                                            precursor_col=parameters_dict["precursor_col"],
                                                            protein_col=parameters_dict["protein_col"],
                                                            on_filter_mask=parameters_dict["define_protein_missingness_on_filter_mask"])
                    #Logging
                    protmissing_counter = Counter(protmissing_mask_df.to_numpy().flatten())
                    log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tDefined missingness if no precursor of a protein was measured as " + parameters_dict["missingness_type_protein_absence"] + "\n"
                    log += "\t\t\t\t" +  str(np.round(protmissing_counter[parameters_dict["missingness_type_protein_absence"]]/sum(protmissing_counter.values())*100,2)) + "% of values are defined as missing due to protein absence.\n"
                helper.logger(log_file, module_output=log)
            else:
                protmissing_mask_df = pd.DataFrame(index=data_df.index,columns=data_df.columns,data="")
            if parameters_dict["eval_missingness_by_detection_threshold"]:
                
                #in pandas version 2.2.2
                #C:\Users\boris\AppData\Local\Temp\ipykernel_18288\907775755.py:233: FutureWarning: Downcasting behavior in `replace` is deprecated and will be removed in a future version. To retain the old behavior, explicitly call `result.infer_objects(copy=False)`. To opt-in to the future behavior, set `pd.set_option('future.no_silent_downcasting', True)`
                #  missingness_mask_df.replace({False: 0, True: 1},inplace=True)
                
                
                #identify hierarchy to evaluate missingness, evaluate missingness, store in one matrix
                if parameters_dict["level_for_intensity_dependent_missingness_estimation"] is not None:
                    levels = metadata_df[parameters_dict["level_for_intensity_dependent_missingness_estimation"]].unique()
                    mar_mnar_mask_df = pd.DataFrame()
                    regression_df = pd.DataFrame()
                    decision_boundary_df = pd.DataFrame(index=levels,columns=["decision_boundary","num_samples"])
                    decision_boundary_df.index.name = parameters_dict["level_for_intensity_dependent_missingness_estimation"]
                    # Logging
                    log =  time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tIntensity dependent missingness estimated per " + parameters_dict["level_for_intensity_dependent_missingness_estimation"] + ". In total " + str(len(levels)) + " levels found in data.\n"
                    helper.logger(log_file, module_output=log)
                    for level in levels:
                        samples_in_level = list(set(data_df.columns) & set(metadata_df[metadata_df[parameters_dict["level_for_intensity_dependent_missingness_estimation"]]== level][parameters_dict["metadata_sample_column"]]))
                        mar_mnar_mask_df_level, regression_df_level, decision_boundary_level = intensity_dependent_missingness_estimation(data_df[samples_in_level],  
                                                                                                                                          axis=parameters_dict["axis"], 
                                                                                                                                          boundary=parameters_dict["intensity_dependent_missingness_decision_boundary"], 
                                                                                                                                          detection_proportion_cap=parameters_dict["intensity_dependent_missingness_detection_proportion_cap"], 
                                                                                                                                          method=parameters_dict["intensity_dependent_missingness_fit_method"], 
                                                                                                                                          n_knots=parameters_dict["intensity_dependent_missingness_n_knots"])
                        #combine mask and regression_df with current level
                        mar_mnar_mask_df = pd.concat([mar_mnar_mask_df, mar_mnar_mask_df_level],axis=1)
                        regression_df_level.columns = [x+"_"+str(level) for x in regression_df_level.columns]
                        regression_df = regression_df.join(regression_df_level,how="outer")
                        decision_boundary_df.loc[level,"decision_boundary"] = decision_boundary_level
                        decision_boundary_df.loc[level,"num_samples"] = len(samples_in_level)
                    # Logging of decision_boundary_dataframe
                    log = "\t\t\t\t" + decision_boundary_df.to_string(index=False).replace("\n","\n\t\t\t\t")+"\n"
                    helper.logger(log_file, log)
                    if parameters_dict["save_intermediate_output_for_QA"]:
                        decision_boundary_df.to_csv(out_dir_fp + "intensity_dependent_missingness_estimation_decision_boundaries.tsv",sep="\t")
                else:
                    # Logging
                    log =  time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tIntensity dependent missingness estimated over the entire dataset.\n"
                    helper.logger(log_file, module_output=log)
                    mar_mnar_mask_df, regression_df, decision_boundary = intensity_dependent_missingness_estimation(data_df, 
                                                                                                                    axis=parameters_dict["axis"], 
                                                                                                                    boundary=parameters_dict["intensity_dependent_missingness_decision_boundary"], 
                                                                                                                    detection_proportion_cap=parameters_dict["intensity_dependent_missingness_detection_proportion_cap"], 
                                                                                                                    method=parameters_dict["intensity_dependent_missingness_fit_method"], 
                                                                                                                    n_knots=parameters_dict["intensity_dependent_missingness_n_knots"])
                    #Logging decision boundary
                    log =  time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tDecision boundary for intensity dependent missingness estimation: " + str(np.round(decision_boundary,2)) + "\n"
                    helper.logger(log_file, module_output=log)
                if parameters_dict["save_intermediate_output_for_QA"]:
                    regression_df.to_csv(out_dir_fp + "intensity_dependent_missingness_estimation_plot_data.tsv",sep="\t")
            else:
                mar_mnar_mask_df = pd.DataFrame(index=data_df.index,columns=data_df.columns,data="")
            # combine masks
            imputation_mask_df = protmissing_mask_df.loc[data_df.index,data_df.columns] + "/" + mar_mnar_mask_df.loc[data_df.index,data_df.columns]
            #replace / in beginning or end
            imputation_mask_df = imputation_mask_df.replace(["^/","/$"],"",regex=True)
            #reset index to precursor id
            imputation_mask_df.index = imputation_mask_df.index.get_level_values(parameters_dict["precursor_col"])
            if parameters_dict["save_intermediate_output_for_QA"]:
                imputation_mask_df.to_csv(out_dir_fp + "imputation_mask.tsv",sep="\t")
                log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tImputation mask storing missingness types created and saved.\n"
            else:
                log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tImputation mask storing missingness types created.\n"
            #Logging missing value types
            log += time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tMissingness types and corresponding numbers of values:\n\t\t\t\t" + str(Counter(imputation_mask_df.to_numpy().flatten())) + "\n"
            helper.logger(log_file, module_output=log)

            #extract missingness types
            if parameters_dict["remember_missingness_type_in_mask"]:
                missingness_types_all = pd.Series(imputation_mask_df.to_numpy().flatten()).replace("",np.nan).dropna().unique()
                missingness_types_df = pd.DataFrame(columns=["filter","value"])
                missingness_types_df["filter"] = missingness_types_all
                missingness_types_df["value"] = [x/10 for x in range(1,(len(missingness_types_all)+1))]
                mask_description_df = pd.concat([mask_description_df,missingness_types_df.set_index("filter")],axis=0)
            
            #remember which and how many values were missing initially
            total_num_values0,num_missing_values0,missing_val_locs0 = obtain_missingness_for_log(data_df=data_df)
            # hierarchical imputation
            if parameters_dict["hierarchical_imputation"]:
                #IMPORTANT: hierarchy values should be unique, e.g. there cannot be plate "1" in several batches, or name them as "batch1_plate1" and "batch2_plate2" etc. by concatenating hierarchies
                sample_hierarchies = impu_table_df.columns.tolist()
                sample_hierarchies = sample_hierarchies[:sample_hierarchies.index("missingness_thresh")] #sample hierarchy columns must be present before column precursor_hierarchy
                sample_hierarchies.reverse() #sort from lowest to highest hierarchy for correct imputation
                for sample_hierarchy in sample_hierarchies: #start with lowest sample hierarchy, go up, e.g. plate, batch, superbatch. This way one can define different precursor_hierarchies such as missingness thresholds and impute on larger hierarchy with higher threshold
                    # select hierarchy values
                    sample_hierarchy_values = impu_table_df[sample_hierarchy].dropna().unique()
                    #Logging
                    log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tImputing per sample hierarchy: " + str(sample_hierarchy) + ".\n"
                    log += "\t\t\t\tThis sample hierarchy contains " + str(len(sample_hierarchy_values)) + " instances.\n"
                    helper.logger(log_file, module_output=log)
                    for sample_hierarchy_value in sample_hierarchy_values:
                        samples = list(set(metadata_df[metadata_df[sample_hierarchy]==sample_hierarchy_value][parameters_dict["metadata_sample_column"]]) & set(data_df.columns))
                        #Logging number of samples
                        log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tImputing per " + str(sample_hierarchy_value) + " on " + str(len(samples)) + " samples.\n"
                        helper.logger(log_file, module_output=log)
                        # select precursor_hierarchy (e.g. missingness 0.3, 0.5)
                        precursor_hierarchies = impu_table_df[impu_table_df[sample_hierarchy]==sample_hierarchy_value]["missingness_thresh"].unique()
                        # calculate missingness per precursor for given samples
                        missingness = data_df[samples].isnull().mean(axis=parameters_dict["axis"]) #if one missingness_type is in imputed in inner loop, it shall not affect missingness calculation. However the precursor could already be imputed.
                        for precursor_hierarchy in precursor_hierarchies:
                            # select precursors matching hierarchy (with prior calculated missingness)
                            precursors = missingness[missingness<precursor_hierarchy].index
                            #Logging number of precursors
                            log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tImputing " + str(len(precursors)) + " precursors with missingness up to " + str(precursor_hierarchy*100) + "%.\n"
                            helper.logger(log_file, module_output=log)
                            # extract missingness types from mask
                            missingness_types = pd.Series(imputation_mask_df.loc[precursors,samples].to_numpy().flatten()).replace("",np.nan).dropna().unique()
                            if len(missingness_types) == 0:
                                log = "\t\t\t\tNo missing values in this sample and precursor subset.\n"
                                helper.logger(log_file, module_output=log)
                                continue
                            for missingness_type in missingness_types:
                                imputation_method = impu_table_df.loc[(impu_table_df[sample_hierarchy]==sample_hierarchy_value) & (impu_table_df["missingness_thresh"]==precursor_hierarchy) & (impu_table_df["missingness_type"]==missingness_type),"imputation_method"].iloc[0]
                                imputation_parameters = impu_table_df.loc[(impu_table_df[sample_hierarchy]==sample_hierarchy_value) & (impu_table_df["missingness_thresh"]==precursor_hierarchy) & (impu_table_df["missingness_type"]==missingness_type)].iloc[0].to_dict()
                                # Impute all missing values for selected precursors and samples, set those which should not be imputed because of different missingness type to nan after imputation
                                # precursors and samples with all missing values in this subset are not imputed at all in this step
                                data_part_df = data_df.loc[precursors,samples].dropna(axis=1,how="all").dropna(axis=0,how="all")
                                #remember which and how many values were missing prior imputation with given missingness_type
                                missing_val_locs_m = getIndexes(imputation_mask_df.loc[data_part_df.index,data_part_df.columns] == missingness_type,True)
                                total_num_values,num_missing_values,missing_val_locs = obtain_missingness_for_log(data_df=data_part_df) # all missingness_types
                                missing_val_locs_m = list(set(missing_val_locs_m) & set(missing_val_locs)) # to exclude prior imputed values
                                num_missing_values_m = len(missing_val_locs_m) 
                                
                                #IMPUTATION FUNCTION CALLS
                                if imputation_method == "knn":
                                    data_part_df = knn_imputation(data_df=data_part_df,k=int(imputation_parameters["k"]),weights=parameters_dict["knn_weights"],metric=parameters_dict["knn_metric"],transposed=True)
                                elif imputation_method == "knn-ML":
                                    data_part_df = knn_imputation(data_df=data_part_df,k=int(imputation_parameters["k"]),weights=parameters_dict["knn_weights"],metric=parameters_dict["knn_metric"],transposed=False)
                                elif imputation_method == "minval-noise":
                                    data_part_df = minval_imputation(data_df=data_part_df,axis=parameters_dict["axis"],add_noise=True)
                                elif imputation_method == "minval":
                                    data_part_df = minval_imputation(data_df=data_part_df,axis=parameters_dict["axis"],add_noise=False)
                                elif imputation_method == "meanval-noise":
                                    data_part_df = meanval_imputation(data_df=data_part_df,axis=parameters_dict["axis"],add_noise=True)
                                elif imputation_method == "meanval":
                                    data_part_df = meanval_imputation(data_df=data_part_df,axis=parameters_dict["axis"],add_noise=False)
                                elif imputation_method == "medianval-noise":
                                    data_part_df = medianval_imputation(data_df=data_part_df,axis=parameters_dict["axis"],add_noise=True)
                                elif imputation_method == "medianval":
                                    data_part_df = medianval_imputation(data_df=data_part_df,axis=parameters_dict["axis"],add_noise=False)
                                elif imputation_method == "bpca":
                                    data_part_df = bpca_fill(data_part_df.T).T
                                elif imputation_method == "graphbased":
                                    data_part_df = graphbased_imputation(data_df=data_part_df,
                                                                         id_mapping_df=id_mapping_df,
                                                                         protein_col=parameters_dict["protein_col"],
                                                                         precursor_col=parameters_dict["precursor_col"],
                                                                         impute_modifications=parameters_dict["graphbased_imputation_impute_modifications"],
                                                                         estimate_threshold=parameters_dict["graphbased_imputation_threshold"],
                                                                         weighting=parameters_dict["graphbased_imputation_weighting"])
    
                                #only allow imputed values for selected missingness type
                                #missingness_mask_df = imputation_mask_df.loc[precursors,samples]
                                #data_part_df = data_part_df[(missingness_mask_df == missingness_type) | (missingness_mask_df == "")]
                                data_part_df =  data_part_df[imputation_mask_df.loc[data_part_df.index,data_part_df.columns] == missingness_type]
                                #update data_df with data_part_df
                                data_df.loc[data_part_df.index,data_part_df.columns] = data_df.loc[data_part_df.index,data_part_df.columns].add(data_part_df,fill_value=0).replace(0,np.nan)
                                # remember which and how many values are missing after imputation
                                still_total_num_values,num_still_missing_values,still_missing_val_locs = obtain_missingness_for_log(data_df=data_part_df)
                                #update mask_df with imputed values
                                if parameters_dict["remember_missingness_type_in_mask"]:
                                    mask_value = mask_description_df.loc[missingness_type,"value"]
                                else:
                                    mask_value = eval("0."+str(int(mask_description_df.loc["Imputed","value"]))) #adjust here if changed Imputed changed in mask_description_df to 0.99 #TODO
                                impu_val_locs = list(set(missing_val_locs_m) - set(still_missing_val_locs))
                                mask_df = update_mask(mask_df=mask_df, locs=impu_val_locs, mask_value=mask_value)
                                #Logging number of imputed values
                                if num_missing_values_m !=0:
                                    frac_imputed_missingness_type = len(impu_val_locs)/num_missing_values_m
                                else:
                                    frac_imputed_missingness_type = 0
                                if num_missing_values!=0:
                                    frac_missingness_type_from_missing = num_missing_values_m/num_missing_values
                                else:
                                    frac_missingness_type_from_missing = 0
                                frac_missingness_type_from_all = num_missing_values_m/total_num_values
                                log = "\t\t\t\t" + str(np.round(frac_missingness_type_from_all*100,2)) + "% of values were missing, of these " + str(np.round(frac_missingness_type_from_missing*100,2)) +"% with missingness type " + missingness_type + "\n"
                                log += "\t\t\t\t" + str(np.round(frac_imputed_missingness_type*100,2)) + "% of values with missingness type " + str(missingness_type) + " were imputed using imputation method " + imputation_method + ".\n" 
                                helper.logger(log_file, module_output=log)

            if not parameters_dict["hierarchical_imputation"] or parameters_dict["impute_over_whole_dataset"]: #no hierarchical imputation
                log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tImputing over entire dataset.\n"
                missingness_types = pd.Series(imputation_mask_df.to_numpy().flatten()).replace("",np.nan).dropna().unique()
                notimplemented_missingness_types = set(missingness_types) - {"MAR","MNAR",parameters_dict["missingness_type_protein_absence"],parameters_dict["missingness_type_protein_absence"] + "/MAR",parameters_dict["missingness_type_protein_absence"] + "/MNAR"}
                if  notimplemented_missingness_types != set():
                    log = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\t###WARNING### No imputation method defined for missingness type(s) " + str(notimplemented_missingness_types) + "\n\t\t\t\tPlease define a variable specifying the missingness_type and change the source code.\n\t\t\t\tMissingness types skipped for now.\n"
                for missingness_type in missingness_types:
                    if missingness_type == "MAR":
                        imputation_method = parameters_dict["imputation_method_MAR"]
                    elif missingness_type == "MNAR":
                        imputation_method = parameters_dict["imputation_method_MNAR"]
                    elif missingness_type == parameters_dict["missingness_type_protein_absence"]:
                        imputation_method = parameters_dict["imputation_method_ProtAbsence"]
                    elif missingness_type == parameters_dict["missingness_type_protein_absence"] + "/MAR":
                        imputation_method = parameters_dict["imputation_method_ProtAbsenceMAR"]
                    elif missingness_type == parameters_dict["missingness_type_protein_absence"] + "/MNAR":
                        imputation_method = parameters_dict["imputation_method_ProtAbsenceMNAR"]
                    else:
                        continue
                    #remember which and how many values were missing prior imputation
                    missing_val_locs_m = getIndexes(imputation_mask_df == missingness_type,True)
                    total_num_values,num_missing_values,missing_val_locs = obtain_missingness_for_log(data_df=data_df) # all missingness_types
                    missing_val_locs_m = list(set(missing_val_locs_m) & set(missing_val_locs))
                    num_missing_values_m = len(missing_val_locs_m)

                    #IMPUTATION FUNCTION CALLS
                    if imputation_method == "knn":
                        data_df = knn_imputation(data_df=data_df,k=parameters_dict["knn_k"],weights=parameters_dict["knn_weights"],metric=parameters_dict["knn_metric"],transposed=True)
                    elif imputation_method == "knn-ML":
                        data_df = knn_imputation(data_df=data_df,k=parameters_dict["knn_k"],weights=parameters_dict["knn_weights"],metric=parameters_dict["knn_metric"],transposed=False)
                    elif imputation_method == "minval-noise":
                        data_df = minval_imputation(data_df=data_df,axis=parameters_dict["axis"],add_noise=True)
                    elif imputation_method == "minval":
                        data_df = minval_imputation(data_df=data_df,axis=parameters_dict["axis"],add_noise=False)
                    elif imputation_method == "meanval-noise":
                        data_df = meanval_imputation(data_df=data_df,axis=parameters_dict["axis"],add_noise=True)
                    elif imputation_method == "meanval":
                        data_df = meanval_imputation(data_df=data_df,axis=parameters_dict["axis"],add_noise=False)
                    elif imputation_method == "medianval-noise":
                        data_df = medianval_imputation(data_df=data_df,axis=parameters_dict["axis"],add_noise=True)
                    elif imputation_method == "medianval":
                        data_df = medianval_imputation(data_df=data_df,axis=parameters_dict["axis"],add_noise=False)
                    elif imputation_method == "bpca":
                        data_df = bpca_fill(data_df.T).T
                    elif imputation_method == "graphbased":
                        data_df = graphbased_imputation(data_df=data_df,
                                                        id_mapping_df=id_mapping_df,
                                                        protein_col=parameters_dict["protein_col"],
                                                        precursor_col=parameters_dict["precursor_col"],
                                                        impute_modifications=parameters_dict["graphbased_imputation_impute_modifications"],
                                                        estimate_threshold=parameters_dict["graphbased_imputation_threshold"],
                                                        weighting=parameters_dict["graphbased_imputation_weighting"])
                    # remember which and how many values are missing after imputation
                    still_total_num_values,num_still_missing_values,still_missing_val_locs = obtain_missingness_for_log(data_df=data_df)
                    #update mask_df with imputed values
                    if parameters_dict["remember_missingness_type_in_mask"]:
                        mask_value = mask_description_df.loc[missingness_type,"value"]
                    else:
                        mask_value = eval("0."+str(mask_description_df.loc["Imputed","value"])) #adjust here if changed Imputed changed in mask_description_df to 0.99 #TODO
                    impu_val_locs = list(set(missing_val_locs_m) - set(still_missing_val_locs))
                    mask_df = update_mask(mask_df=mask_df, locs=impu_val_locs, mask_value=mask_value)
                    #Logging number of imputed values
                    if num_missing_values_m != 0:
                        frac_imputed_missingness_type = len(impu_val_locs)/num_missing_values_m
                    else:
                        frac_imputed_missingness_type = 0
                    if num_missing_values != 0:
                        frac_missingness_type_from_missing = num_missing_values_m/num_missing_values
                    else:
                        frac_missingness_type_from_missing = 0
                    frac_missingness_type_from_all = num_missing_values_m/total_num_values
                    log = "\t\t\t\t" + str(np.round(frac_missingness_type_from_all*100,2)) + "% of values were missing, of these " + str(np.round(frac_missingness_type_from_missing*100,2)) +"% with missingness type " + missingness_type + "\n"
                    log += "\t\t\t\t" + str(np.round(frac_imputed_missingness_type*100,2)) + "% of values with missingness type " + str(missingness_type) + " were imputed using imputation method " + imputation_method + ".\n" 

                    helper.logger(log_file, log)
            #Count final number of imputed and missing values
            num_values,num_still_missing_values,still_missing_val_locs = obtain_missingness_for_log(data_df=data_df)
            frac_imputed_values = (num_missing_values0 - num_still_missing_values) / total_num_values0
            log = time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputation finished.\n"
            log += "\t\t\t\tFinally " + str(np.round(frac_imputed_values*100,2)) + "% of the dataset was imputed.\n"
            log += "\t\t\t\t" + str(np.round((num_missing_values0-num_still_missing_values)/num_missing_values0*100,2)) + "% of the missing values were imputed.\n"
            log += "\t\t\t\t" + str(np.round(num_still_missing_values/total_num_values0*100,2)) + "% of the dataset remains missing.\n"
            helper.logger(log_file,log)
            

        except Exception as e:
            helper.logger(log_file, "###ERROR### Module ran into an error: " + e)
            sys.exit("Imputation exited with an error.")
    else:
        log = time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tNo missing values identified, imputation skipped.\n"
        helper.logger(log_file,log)
    
    module_output = time.strftime("%Y-%m-%d%t%H:%M:%S") + "\tmodule_imputation ends without errors.\n"  
    helper.close_module(data_processed_df=data_df, masked_df=mask_df, mask_description_df=mask_description_df, log_file=log_file, module_output=module_output, out_dir_fp=out_dir_fp, outname=outname, mask=True)
       
    
