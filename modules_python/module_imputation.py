"""
Input: wide-format data
Output: imputed wide-format data

    - impute_missing_peptides: self-made graph-based imputation strategy
    - getIndexes: returns a list of positions where element exists in the given dataframe
    - knn_imputation: performs k-nearest-neighbour imputation on dataframe
    - minval_imputation: performs modified minimum value imputation on dataframe
    - meanval_imputation: performs mean value imputation on dataframe
    - bpca_imputation: performs Bayesian PCA imputation on dataframe



"""

#general
import pandas as pd
import numpy as np
import re
import sys
import time

#own scripts
from Python.Snippets import global_variables as global_variables # left this as placeholder for logging, when adding logging to main functions this will be moved to main, can be commented
from modules_python import module_helper

#statistics and modelling
from itertools import permutations
import igraph
from sklearn.impute import KNNImputer
import statsmodels.api as sm

##############################################################################
# Mixed imputation by detection probability modelling using precursor intensities
# Author: SB
##############################################################################

def preprocess_data(data):
    """
    Preprocesses input data by removing fully missing columns and replacing zeros with NaNs.

    author: SB

    Parameters
    ----------
    data : pd.DataFrame
        Input data for preprocessing.

    Returns
    -------
    pd.DataFrame
        Preprocessed DataFrame.
    """

    data = data.dropna(axis=1, how='all')
    data = data.replace(0, np.nan)

    return data

def compute_log2_means_and_missingness(data):
    """
    Computes log2-transformed mean abundances and calculates missingness proportion per precursor.

    author: SB

    Parameters
    ----------
    data : pd.DataFrame
        Input DataFrame.

    Returns
    -------
    means_log2 : pd.Series
        Log2-transformed mean abundance per precursor.
    missing : pd.Series
        Missingness proportion per precursor.
    """

    means = data.mean(axis=0)
    means_log2 = np.log2(means.replace(0, np.nan)).replace([np.inf, -np.inf], np.nan)
    missing = data.isnull().mean(axis=0)
    return means_log2, missing

def detection_probability_curve(means_log2, missingness, boundary=0.5):
    """
    Fits logistic regression to precursor detection probability and calculates a decision 
    boundary for missingness classification.

    author: SB

    Parameters
    ----------
    means_log2 : pd.Series
        Log2-transformed precursor abundances.
    missingness : pd.Series
        Missingness proportion per precursor.
    boundary : float, optional
        Decision probability boundary for classification (default 0.5).

    Returns
    -------
    decision_boundary : float
        Calculated abundance threshold for MNAR/MAR classification.
    result : statsmodels.GLMResults
        Logistic regression model fit results.
    """

    df = pd.DataFrame({
        'Avg Log2 Abundance': means_log2,
        'Missing %': 1 - missingness
    }).dropna()

    X = sm.add_constant(df['Avg Log2 Abundance'])
    y = df['Missing %']

    model = sm.GLM(y, X, family=sm.families.Binomial())
    result = model.fit()

    intercept, coef = result.params

    decision_boundary = (np.log(boundary / (1 - boundary)) - intercept) / coef

    return round(decision_boundary, 2), result


def classify_precursors(means_log2, boundary):
    """
    Classifies precursors as Missing Not At Random (MNAR) based on an abundance boundary.

    author: SB

    Parameters
    ----------
    means_log2 : pd.Series
        Log2 mean abundances.
    boundary : float
        Decision boundary for classification.

    Returns
    -------
    pd.Series
        Boolean mask indicating MNAR (True) or MAR (False) precursors.
    """

    MNAR_mask = means_log2 < boundary
    return MNAR_mask


def mixed_imputer(data, MNAR_mask, missing, boundary, ramp, knn=3):
    """
    Imputes missing data using minimum abundance for MNAR precursors and KNN imputation for MAR precursors.

    author: SB

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame with missing precursor abundances.
    MNAR_mask : pd.Series
        Boolean mask identifying MNAR precursors.
    missing : pd.Series
        Missingness per precursor.
    boundary : float
        Decision boundary for MNAR/MAR classification.
    ramp : float
        Missingness threshold percentage for precursor inclusion.
    knn : int, optional
        Number of neighbors for KNN imputation (default 3).

    Returns
    -------
    pd.DataFrame
        DataFrame with missing values imputed.
    """

    threshold = ramp / 100
    min_values = data.min()

    # Filter columns based on missingness threshold
    data_filtered = data.loc[:, missing < (1 - threshold)]

    # Ensure MNAR and MAR columns exist after filtering
    MNAR_cols = MNAR_mask[MNAR_mask].index.intersection(data_filtered.columns)
    MAR_cols = MNAR_mask[~MNAR_mask].index.intersection(data_filtered.columns)

    # Minimum imputation for MNAR columns
    data_filtered[MNAR_cols] = data_filtered[MNAR_cols].fillna(min_values[MNAR_cols])

    # KNN imputation (applied after minimum imputation)
    imputer = KNNImputer(n_neighbors=knn)
    data_imputed = pd.DataFrame(
        imputer.fit_transform(data_filtered),
        index=data_filtered.index,
        columns=data_filtered.columns
    )

    return data_imputed

def mixed_imputation(data, threshold = 50):
    """
    Integrated pipeline for preprocessing data, classifying precursors, and performing 
    mixed imputation (minimum and KNN).

    author: SB

    Parameters
    ----------
    data : pd.DataFrame
        Input data requiring imputation.
    threshold : float, optional
        Missingness threshold for precursor retention (default 50%).

    Returns
    -------
    pd.DataFrame
        Imputed dataset.
    """

    # Initial Data
    sample_filtered = preprocess_data(data)

    # Compute Means and Missingness
    means_log2, missing = compute_log2_means_and_missingness(sample_filtered)

    # Fit Logistic Regression
    boundary, regression_result = detection_probability_curve(means_log2, missing, boundary=0.5)
    print(f"Decision boundary (Avg Log2): {boundary}")

    # Classify MNAR/MAR
    MNAR_mask = classify_precursors(means_log2, boundary)

    # Perform Mixed Imputation
    imputed_data = mixed_imputer(sample_filtered, MNAR_mask, missing, boundary, ramp=threshold, knn=3)

    return imputed_data

def graphbased_imputation(peptide_profiles,mask=None,mask_description=None,impute_modifications=False,group=None,estimate_threshold=1,weighting=False,keep_track=True):
    """
    Graph-based approach to impute missing precursors. This approach assumes that the ratio of any two precursors from a given protein solely relies on physico-chemical properties and digestion and therefore is equal in all samples/runs.
    Creates for every protein group a directed, weighted graph with precursors as nodes and median ratio of any two precursors (nodes) over all samples/fraction as edge weights between nodes.
    These ratios are then multiplied with present precursors over samples, median over them is used as intensity of missing precursor. 
    This approach may not impute all missing precursors, as advantage it uses global information over all samples about precursor ratios, but imputes missing values only from precursor intensities within a sample times global ratios.
    Imputed values are marked in mask as a second number after decimal point in the matrix, preserving the information if the precursor was previously filtered or not measured at all
    
    author: BZ
    
    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None.
    impute_modifications : bool, optional
        Switch to True to also impute precursors with modifications. In fractionation experiments the ratio between modified precursors is not necessarily equal over fractions(=samples). The default is False.
    group : str or int, optional
        Specifies sample group of given peptide_profiles; it is used to find the correct mapping in the filtering mask. The default is None.
    estimate_threshold : int, optional
        Minimum number of annotated edges (known ratios between two precursors) to impute a missing precursor value. The default is 1.
    weighting : bool, optional
        If True calculates the weighted median between estimated imputed values for a precursor weighted on how often the ratio between two precursors was used to calculate edge weights. The default is False.
    keep_track : bool, optional
        If True, imputed precursors are remembered in mask. The default is True.
        
    Returns
    -------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities with imputed values.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.

    """
    #num_measured_peptides=len(peptide_profiles.reshape(peptide_profiles.shape[0]*peptide_profiles.shape[1],1).dropna())
    num_peptides=peptide_profiles.shape[0]*peptide_profiles.shape[1]
    #num_missing_peptides=np.sum(pd.isna(peptide_profiles.reshape(peptide_profiles.shape[0]*peptide_profiles.shape[1],1)))
    #peptide_profiles.replace(0,np.nan,inplace=True)
    num_missing_peptides=pd.isna(peptide_profiles).sum().sum()
    print("Begin Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputing missing precursors by own graph-based method...\n")
        logfile.close()
    '''
    if "sample group" in mask.index:
        indices = mask.index.drop("sample group")
    else:
        indices = mask.index
    '''
    indices=peptide_profiles.index
    #for every protein group
    for pg in indices.get_level_values("Protein.Group").unique():
        #create graph for protein with precursors as nodes and ratios as edges
        temp_df = peptide_profiles[peptide_profiles.index.get_level_values("Protein.Group")==pg]
        #exclude columns without any precursor intensities for given protein group
        temp_df.dropna(thresh=1,axis=1,inplace=True)
        #exclude precursors with modification sites
        if impute_modifications==False:
            modi_precs={x if "UniMod:" in x else np.nan for x in temp_df.index.get_level_values("Precursor.Id")}
            if np.nan in modi_precs:
                modi_precs.remove(np.nan)
            if modi_precs!=set():
                #include non-modified sequence
                modi_precs=modi_precs | {re.sub("\(.*\)","",prec) for prec in modi_precs}
                temp_df.drop(temp_df[temp_df.index.get_level_values("Precursor.Id").isin(modi_precs)].index,inplace=True)
        if temp_df.empty:
            continue
        #edges=permutations(range(len(temp_df.index)),2)
        #edge_weights=list(map(lambda x: np.nanmedian(temp_df.loc[x[0]]/temp_df.loc[x[1]]),permutations(temp_df.index,2)))
        #ratio_weights tell how often the ratio between two given precursors was used to calculate edge weight. Edges with higher ratio_weight have higher impact in calculation of median for the missing value
        g=igraph.Graph(directed=True)
        g.add_vertices(range(len(temp_df.index)))
        g.vs["name"]=temp_df.index
        g.add_edges(permutations(range(len(temp_df.index)),2))
        #g.es["weight"]=list(map(lambda x: np.nanmedian(temp_df.loc[x[1]].astype(float)/temp_df.loc[x[0]].astype(float)),permutations(temp_df.index,2)))
        g.es["weight"]=list(map(lambda x: np.nanmedian(temp_df.loc[x[1]]/temp_df.loc[x[0]]),permutations(temp_df.index,2)))
        if weighting:
            g.es["ratio_weight"]=list(map(lambda x: len((temp_df.loc[x[1]]/temp_df.loc[x[0]]).dropna()),permutations(temp_df.index,2)))
        g.es["name"]=list(permutations(temp_df.index,2))
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
                        peptide_profiles.loc[precursor,col]=module_helper.weighted_nanmedian(estimates,ratio_weights)
                    else:
                        peptide_profiles.loc[precursor,col]=np.nanmedian(estimates)
                #[g.es.find(name=(x["name"],precursor))["weight"]*temp_df.loc[x["name"],col] for x in g.vs.find(name=precursor).predecessors()]
                #[g.es.find(name=(x["name"],precursor))["weight"] for x in g.vs.find(name=precursor).predecessors()]
                #[temp_df.loc[x["name"],col] for x in g.vs.find(name=precursor).predecessors()]
                #remember imputed precursors in mask
                if keep_track:
                    if peptide_profiles.loc[precursor,col]!=0 & pd.notna(peptide_profiles.loc[precursor,col]):
                        if type(mask.columns)==pd.core.indexes.multi.MultiIndex:
                            mask.loc[precursor,(col,group)]=mask.loc[precursor,(col,group)]+float("0."+str(mask_description.loc["Imputed","value"]))
                        else:
                            mask.loc[precursor,col]=mask.loc[precursor,col]+eval("0."+str(mask_description.loc["Imputed","value"]))
    print("End Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    num_missing_peptides_after_imputation=pd.isna(peptide_profiles).sum().sum()
    print(str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation")
    print(str(np.round(num_missing_peptides_after_imputation/num_peptides*100,2))+"% of precursors are missing after imputation")
    #(num peptides - nans after imputation) / (num_peptides - nans before imputation)
    print(str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed")
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputation finished.\n\t\t\t\t"+str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation"+"\n")
        logfile.write("\t\t\t\t"+str(np.round(num_missing_peptides_after_imputation/num_peptides*100,2))+"% of precursors are missing after imputation\n")
        logfile.write("\t\t\t\t"+str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed ("+str(num_peptides-num_missing_peptides)+" Precursors)\n")
        logfile.close()
    return([peptide_profiles, mask, mask_description])


#copied from https://www.geeksforgeeks.org/find-location-of-an-element-in-pandas-dataframe-in-python/ , access 10.02.2023, 11:44
def getIndexes(dfObj, value):
    """
    This function will return a list of positions where element exists in the dataframe.
    
    author: BZ
    
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


def knn_imputation(peptide_profiles,mask=None,mask_description=None,k=5,group=None,keep_track=True,transposed=False):
    """
    k-nearest neighbours imputation on peptide or protein profiles.

    author: BZ
    
    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None.
    k : int, optional
        number of neighbours to take into account for kNN Imputation. The default is 5.
    group : str or int, optional
        Specifies sample group of given peptide_profiles; it is used to find the correct mapping in the filtering mask. The default is None.
    keep_track : bool, optional
        If True, imputed precursors are remembered in mask. The default is True.
    transposed : bool, optional
        If True, tranposes the data matrix prior to imputation and impute as described by Troyanskaja et al.https://doi.org/10.1093/bioinformatics/17.6.520 , else imputes missing values of a protein based on values of k neighbouring samples


    Returns
    -------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.

    """
    num_peptides=peptide_profiles.shape[0]*peptide_profiles.shape[1]
    #peptide_profiles.replace(0,np.nan,inplace=True)
    num_missing_peptides=pd.isna(peptide_profiles).sum().sum()
    print("Begin Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputing missing precursors using kNN (k="+str(k)+")...\n")
        logfile.close()
    if keep_track:
        #remember which values were missing in mask
        if mask is not None:
            missing_val_locs=getIndexes(peptide_profiles.isna(),True)
    #knn imputation
    imputer = KNNImputer(n_neighbors=k,weights="distance",keep_empty_features=True, metric="nan_euclidean")
    indices=peptide_profiles.index
    columns=peptide_profiles.columns
    if transposed:
        imputed_values = imputer.fit_transform(peptide_profiles)
        peptide_profiles = pd.DataFrame(index=indices,columns=columns,data=imputed_values)
    else:
        imputed_values = imputer.fit_transform(peptide_profiles.T)
        peptide_profiles = pd.DataFrame(index=columns,columns=indices,data=imputed_values).T
    if keep_track:
        #assign imputed precursors in mask
        if type(mask.columns)==pd.core.indexes.multi.MultiIndex and not group in ["",None]:
            for precursor,column in missing_val_locs:
                mask.loc[precursor,(column,group)]=mask.loc[precursor,(column,group)]+float("0."+str(mask_description.loc["Imputed","value"]))
        else:
            for precursor,column in missing_val_locs:
                mask.loc[precursor,column]=mask.loc[precursor,column]+eval("0."+str(mask_description.loc["Imputed","value"]))
    #print and store information in logfile
    print("End Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    num_missing_peptides_after_imputation=pd.isna(peptide_profiles).sum().sum()
    print(str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation")
    #(num peptides - nans after imputation) / (num_peptides - nans before imputation)
    print(str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed")
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputation finished.\n"+str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation"+"\n")
        logfile.write("\t\t\t\t"+str(np.round(num_missing_peptides_after_imputation/num_peptides*100,2))+"% of precursors are missing after imputation\n")
        logfile.write("\t\t\t\t"+str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed ("+str(num_peptides-num_missing_peptides)+" Precursors)\n")
        logfile.close()
    return([peptide_profiles, mask, mask_description])
   

def minval_imputation(peptide_profiles,mask=None,mask_description=None,group=None,keep_track=True):
    """
    minimum value imputation on peptide or protein profiles. Each missing value is imputed by the minimum of the row divided by the square root of 2 times a random number between 0.75 and 1.25.
    
    author: BZ
    
    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None.
    group : str or int, optional
        Specifies sample group of given peptide_profiles; it is used to find the correct mapping in the filtering mask. The default is None.
    keep_track : bool, optional
        If True, imputed precursors are remembered in mask. The default is True.


    Returns
    -------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. 
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.

    """
    num_peptides=peptide_profiles.shape[0]*peptide_profiles.shape[1]
    #peptide_profiles.replace(0,np.nan,inplace=True)
    num_missing_peptides=pd.isna(peptide_profiles).sum().sum()
    print("Begin Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputing missing precursors using minimum values per precursor divided by sqrt(2) times random number [0.75,1.25]...\n")
        logfile.close()
    if keep_track:
        #remember which values were missing in mask
        if mask is not None:
            missing_val_locs=getIndexes(peptide_profiles.isna(),True)
    #minimum value imputation
    #simple imputation
    ##peptide_profiles.replace(np.nan,peptide_profiles.min().min())
    #impute by minimum value per ID divided by sqrt(2) times random number between 0.75 and 1.25
    for row in peptide_profiles.index:
        peptide_profiles.loc[row]=peptide_profiles.loc[row].apply(lambda x: peptide_profiles.loc[row].min()/np.sqrt(2)*np.random.uniform(0.75,1.25) if pd.isna(x) else x)
    if keep_track:
        #assign imputed precursors in mask
        if mask is not None:
            if type(mask.columns)==pd.core.indexes.multi.MultiIndex and group not in ["",None]:
                for precursor,column in missing_val_locs:
                    mask.loc[precursor,(column,group)]=mask.loc[precursor,(column,group)]+float("0."+str(mask_description.loc["Imputed","value"]))
            else:
                for precursor,column in missing_val_locs:
                    mask.loc[precursor,column]=mask.loc[precursor,column]+eval("0."+str(mask_description.loc["Imputed","value"]))
    #print and store information in logfile
    print("End Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    num_missing_peptides_after_imputation=pd.isna(peptide_profiles).sum().sum()
    print(str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation")
    #(num peptides - nans after imputation) / (num_peptides - nans before imputation)
    print(str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed")
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputation finished.\n"+str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation"+"\n")
        logfile.write("\t\t\t\t"+str(np.round(num_missing_peptides_after_imputation/num_peptides*100,2))+"% of precursors are missing after imputation\n")
        logfile.write("\t\t\t\t"+str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed ("+str(num_peptides-num_missing_peptides)+" Precursors)\n")
        logfile.close()
    return([peptide_profiles, mask, mask_description])
    

def meanval_imputation(peptide_profiles,mask=None,mask_description=None,group=None,keep_track=True):
    """
    meam value imputation on peptide or protein profiles. Each missing value is imputed by the mean of the row.
    
    author: BZ
    
    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None.
    group : str or int, optional
        Specifies sample group of given peptide_profiles; it is used to find the correct mapping in the filtering mask. The default is None.
    keep_track : bool, optional
        If True, imputed precursors are remembered in mask. The default is True.


    Returns
    -------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. 
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.

    """
    num_peptides=peptide_profiles.shape[0]*peptide_profiles.shape[1]
    #peptide_profiles.replace(0,np.nan,inplace=True)
    num_missing_peptides=pd.isna(peptide_profiles).sum().sum()
    print("Begin Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputing missing precursors using mean values per precursor...\n")
        logfile.close()
    if keep_track:
        #remember which values were missing in mask
        if mask is not None:
            missing_val_locs=getIndexes(peptide_profiles.isna(),True)
    #minimum value imputation
    #simple imputation
    ##peptide_profiles.replace(np.nan,peptide_profiles.min().min())
    #impute by mean value per ID
    for row in peptide_profiles.index:
        peptide_profiles.loc[row]=peptide_profiles.loc[row].apply(lambda x: peptide_profiles.loc[row].mean() if pd.isna(x) else x)
    if keep_track:
        #assign imputed precursors in mask
        if mask is not None:
            if type(mask.columns)==pd.core.indexes.multi.MultiIndex and group not in ["",None]:
                for precursor,column in missing_val_locs:
                    mask.loc[precursor,(column,group)]=mask.loc[precursor,(column,group)]+float("0."+str(mask_description.loc["Imputed","value"]))
            else:
                for precursor,column in missing_val_locs:
                    mask.loc[precursor,column]=mask.loc[precursor,column]+eval("0."+str(mask_description.loc["Imputed","value"]))
    #print and store information in logfile
    print("End Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    num_missing_peptides_after_imputation=pd.isna(peptide_profiles).sum().sum()
    print(str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation")
    #(num peptides - nans after imputation) / (num_peptides - nans before imputation)
    print(str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed")
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputation finished.\n"+str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation"+"\n")
        logfile.write("\t\t\t\t"+str(np.round(num_missing_peptides_after_imputation/num_peptides*100,2))+"% of precursors are missing after imputation\n")
        logfile.write("\t\t\t\t"+str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed ("+str(num_peptides-num_missing_peptides)+" Precursors)\n")
        logfile.close()
    return([peptide_profiles, mask, mask_description])



#downloaded from https://github.com/unifey/bPCAfill, 18.03.2024 17:08, credit to unifey
def table(x): 
    try: return pd.DataFrame(x)
    except: return pd.DataFrame(list(x.items()))


#downloaded from https://github.com/unifey/bPCAfill, 18.03.2024 17:08, credit to unifey
def bpca_fill(dX, max_iter = 200):
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


def bpca_imputation(peptide_profiles,mask=None,mask_description=None,group=None,keep_track=True):
    """
    Bayesian PCA imputation on peptide or protein profiles as described by Oba et al. https://doi.org/10.1093/bioinformatics/btg287
    
    author: BZ
    
    Parameters
    ----------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. The default is None.
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters. The default is None.
    group : str or int, optional
        Specifies sample group of given peptide_profiles; it is used to find the correct mapping in the filtering mask. The default is None.
    keep_track : bool, optional
        If True, imputed precursors are remembered in mask. The default is True.
    
    
    Returns
    -------
    peptide_profiles : pd.DataFrame
        rows equal precursor ids, columns correspond to samples and values in the dataframe to intensities.
    mask : pd.DataFrame.
        dataframe with precursors as rows and runs as columns, with values explaining why precursors were filtered. 
    mask_description: pd.DataFrame.
        mapping numbers in mask to filters.
    
    """
    num_peptides=peptide_profiles.shape[0]*peptide_profiles.shape[1]
    #peptide_profiles.replace(0,np.nan,inplace=True)
    num_missing_peptides=pd.isna(peptide_profiles).sum().sum()
    print("Begin Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputing missing precursors using BPCA...\n")
        logfile.close()
    if keep_track:
        #remember which values were missing in mask
        if mask is not None:
            missing_val_locs=getIndexes(peptide_profiles.isna(),True)
    #BPCA Imputation
    peptide_profiles=bpca_fill(peptide_profiles.T).T
    if keep_track:
        #assign imputed precursors in mask
        if mask is not None:
            if type(mask.columns)==pd.core.indexes.multi.MultiIndex and group not in ["",None]:
                for precursor,column in missing_val_locs:
                    mask.loc[precursor,(column,group)]=mask.loc[precursor,(column,group)]+float("0."+str(mask_description.loc["Imputed","value"]))
            else:
                for precursor,column in missing_val_locs:
                    mask.loc[precursor,column]=mask.loc[precursor,column]+eval("0."+str(mask_description.loc["Imputed","value"]))
    #print and store information in logfile
    print("End Data Imputation: "+time.strftime("%Y-%m-%d%t%H:%M:%S"))
    num_missing_peptides_after_imputation=pd.isna(peptide_profiles).sum().sum()
    print(str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation")
    #(num peptides - nans after imputation) / (num_peptides - nans before imputation)
    print(str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed")
    if global_variables.logfilename is not None:
        logfile=open(global_variables.logfilename,"a")
        logfile.write(time.strftime("%Y-%m-%d%t%H:%M:%S")+"\tImputation finished.\n"+str(np.round(num_missing_peptides/num_peptides*100,2))+"% of precursors were missing before imputation"+"\n")
        logfile.write("\t\t\t\t"+str(np.round(num_missing_peptides_after_imputation/num_peptides*100,2))+"% of precursors are missing after imputation\n")
        logfile.write("\t\t\t\t"+str(np.round(100-(num_peptides-num_missing_peptides)/(num_peptides-num_missing_peptides_after_imputation)*100,2))+"% of precursors are imputed ("+str(num_peptides-num_missing_peptides)+" Precursors)\n")
        logfile.close()
    return([peptide_profiles, mask, mask_description])
