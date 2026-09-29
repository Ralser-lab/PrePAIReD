#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 16:00:21 2025

@author: oliverlemke
"""

import sys
import os
import argparse
import re
import ast
from tkinter.filedialog import askopenfilename, askdirectory

import numpy as np
import pandas as pd
import statsmodels.api as sm

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter

from sklearn.decomposition import PCA
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.preprocessing import StandardScaler, SplineTransformer

from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr, pointbiserialr, chi2

import module_helper as helper


### Generalize colors as an input_variable
### PCA: Option to plot centroids
### Scatter: Dot size as input parameter

### Possible extensions
## Options to add later
## Group-wise plots
## if sample list is provided samples vs. all
## Clustermap of DIANN metrics
## Clustermap with groups (before and after), multi layer (QC vs. biological)

## More comments for main

def log_legend(captions, key_figure, figure_path, addition=""):
    """
    Automatic logging of the captions for each figre

    Parameters
    ----------
    captions : str
        Default caption for figure with name key_figure.
    key_figure : str
        Descriptor for the figure. Used to automatically call the corresponding caption.
    figure_path : str
        Path were the figure key_figure will be saved.
    addition : str, optional
        Additional description like selected columns or coloring. The default is "".
    """
    
    log = open(os.path.join(re.search(f"(.+?){key_figure}",figure_path).group(1),"Figure_log.txt"),"a")
    log.write("###FIGURE###: "+figure_path+"\n")
    log.write("###CAPTION###: "+captions.loc[key_figure]["Description"]+addition+"\n")
    log.close()

def get_CV_log(values, base=2):
    """
    Calculate the coefficient of variation (CV) assuming a log-normal distribution and logged data.

    Parameters
    ----------
    values : np.ndarray 
        Array of values for which the CV should be calculated.
    base : float, optional
        Base of the logarithmic transformation. The default is 2.

    Returns
    -------
    cv : float
        Coefficient of variation for provided values.
    """
    
    cv = np.sqrt(float(base**(np.std(values, ddof=1)**2)-1))
    return cv

def get_CV_log_array(df, base=2):
    """
    Calculate the coefficient of variation (CV) assuming a log-normal distribution and logged data for whole data frames (precursors x samples).

    Parameters
    ----------
    values : np.ndarray 
        Array of values for which the CV should be calculated.
    base : float, optional
        Base of the logarithmic transformation. The default is 2.

    Returns
    -------
    cv : list
        Coefficient of variation for provided values.
    """
    
    cv = np.sqrt(np.float64(base**(np.std(df, ddof=1, axis=1)**2)-1))
    return list(cv)

def get_indices_by_type(data_origin, metadata, data_processed=None, sample_type_column=None, sample_type = None):
    """
    Get indices of selected samples according to a sample type e.g. MS.QC-Samples.

    Parameters
    ----------
    data_origin : pd.DataFrame
        Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
    metadata : pd.DataFrame
        Metadata for each sample with the same indices as the columns of data_origin.
    data_processed : pd.DataFrame, optional
        Data (Indices: Precursors, Columns: Samples) after the preprocessing step.
    sample_type_column : str, optional
        Column in Metadata that distinguishes sample type. If None, samples will not be distinguished. The default is None.
    sample_type : str, optional
        Sample type in Metadata[sample_type_column] for which indices should be extracted. The default is None.

    Returns
    -------
    indices_origin : list of str
        Indices of selected sample in the data before the preprocessing step.
    indices_processed: list of str, optional
        Indices of selected sample in the data after the  preprocessing step.
    """

    if sample_type_column:
        if not sample_type:
            raise ValueError("sample_type is not defined.")
        indices_origin = [index for index in metadata.index if metadata.loc[index,sample_type_column]==sample_type and index in data_origin.columns]
        if data_processed is not None:
            indices_processed = [index for index in metadata.index if metadata.loc[index,sample_type_column]==sample_type and index in data_processed.columns]
    else:
        indices_origin = [index for index in metadata.index if index in data_origin.columns]
        if data_processed is not None:
            indices_processed = [index for index in metadata.index if index in data_processed.columns]
    
    if not indices_origin:
        raise ValueError("No indices for selection in original found.")
    if data_processed is not None:
        if not indices_processed:
            raise ValueError("No indices for selection in processed found.")
    
    if data_processed is not None:
        return indices_origin, indices_processed
    else:
        return indices_origin

def get_colors_by_metadata(metadata, column_color, indices_origin, indices_processed=None, column_values_origin=None, column_values_processed=None, cmap="viridis", continuous=False):
    """
    

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for each sample with the same indices as the columns of data_origin.
    column_color : str
        Column in Metadata that is used for coloring the data.
    indices_origin : list of str
        Indices of selected sample in the data before the preprocessing step.
    indices_processed: list of str, optional
        Indices of selected sample in the data after the preprocessing step.
    column_values_origin : pd.Series, optional
        Color values per sample for indices_origin. If set to None column values are selected according to metadata.loc[indices_origin,column_color]. The default is None.
    column_values_processed : pd.Series, optional
        Color values per sample for indices_processed. If set to None column values are selected according to metadata.loc[indices_processed,column_color]. The default is None.
    cmap : str, optional
        Colormap used for coloring. The default is "viridis".
    continuous : bool, optional
        If True continuous data are assumed for coloring e.g. Completeness. If False discrete data are assumed e.g. Batch. The default is False.

    Returns
    -------
    colors_origin : list of tuple of float
        Color values RGBA per sample before preprocessing step.
    colors_processed : list of tuple of float, optional
        Color values RGBA per sample before preprocessing step.
    norm : matplotlib.colors.Normalize
        Norm for the color scale
    values_dict : dict
        Dictionary for mapping values to color. Only filled if discrete values are used for coloring.
    """
    
    if column_values_origin is None:
        if column_color not in metadata.columns:
            raise ValueError("Column not found in metadata")
        column_values_origin = metadata.loc[indices_origin,column_color]
        if indices_processed is not None:
            column_values_processed = metadata.loc[indices_processed,column_color]
            
    if indices_processed is not None:            
        colors_origin, colors_processed, norm, values_dict = get_colors_selection(column_values_origin, column_values_processed, cmap=cmap, continuous=continuous)
        return colors_origin, colors_processed, norm, values_dict
    else:
        colors_origin, norm, values_dict = get_colors_selection(column_values_origin, column_values_processed=None, cmap=cmap, continuous=continuous)
        return colors_origin, norm, values_dict

def get_colors_selection(column_values_origin, column_values_processed=None, cmap="viridis", continuous=False):
    """
    Get color values based for a continuous or discrete pd.Series of entries.

    column_values_origin : pd.Series, optional
        Property value per sample before the preprocsessing step.
    column_values_processed : pd.Series, optional
        Property value per sample after the preprocsessing step.    
    cmap : str, optional
        Colormap used for coloring. The default is "viridis".
    continuous : bool, optional
        If True continuous data are assumed for coloring e.g. Completeness. If False discrete data are assumed e.g. Batch. The default is False.

    Returns
    -------
    colors_origin : list of tuple of float
        Color values RGBA per sample before preprocessing step.
    colors_processed : list of tuple of float, optional
        Color values RGBA per sample before preprocessing step.
    norm : matplotlib.colors.Normalize
        Norm for the color scale
    values_dict : dict
        Dictionary for mapping values to color. Only filled if discrete values are used for coloring.
    """

    cmap_ = plt.get_cmap(cmap)
    
    if continuous:
        if column_values_processed is not None:
            vmin = min(min(column_values_origin),min(column_values_processed))
            vmax = max(max(column_values_origin),max(column_values_processed))
        else:
            vmin = min(column_values_origin)
            vmax = max(column_values_origin)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        values_dict = {}
    else:
        indexes = np.unique(column_values_origin, return_index=True)[1]
        if isinstance(column_values_origin, pd.Series):
            values = [column_values_origin.values[index] for index in sorted(indexes)]
        else:
            values = [column_values_origin[index] for index in sorted(indexes)]
        values_dict = {value:(1+2*ind)/(2*len(values)) for ind,value in enumerate(values)}
        column_values_origin = np.vectorize(values_dict.get)(column_values_origin)
        bounds = np.linspace(0,1,len(values)+1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap_.N)
        if column_values_processed is not None:
            column_values_processed = np.vectorize(values_dict.get)(column_values_processed)
        
    colors_origin = [cmap_(norm(el)) for el in column_values_origin]
    if column_values_processed is not None:
        colors_processed = [cmap_(norm(el)) for el in column_values_processed]
        return colors_origin, colors_processed, norm, values_dict
    return colors_origin, norm, values_dict

def get_color_values_by_property(data_origin, coloring_property, data_processed=None, base=2):
    """
    Calculates the selected property per sample.

    Parameters
    ----------
    data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
    coloring_property : str
        Property used for coloring. Options: "Missingness", Mean_Intensities", "Median_Intensities", Sum_Intensities"
    data_processed : pd.DataFrame, optional
	    Data (Indices: Precursors, Columns: Samples) after the preprocessing step.
    base : float, optional
        Base of the logarithmic transformation. The default is 2.

    Returns
    -------
    column_values_origin : np.ndarray
        Calculated property before the preprocessing step
    column_values_processed : np.ndarray
        Calculated property after the preprocessing step
    """
    
    if coloring_property == "Missingness":
        column_values_origin = np.sum(~np.isnan(data_origin.values),axis=0)
    elif coloring_property=="Mean_Intensities":
        column_values_origin = np.nanmean(data_origin.values,axis=0)
    elif coloring_property=="Median_Intensities":
        column_values_origin = np.nanmedian(data_origin.values,axis=0)
    elif coloring_property=="Sum_Intensities":
        column_values_origin = np.log(np.nansum(base**(data_origin.values),axis=0))/np.log(base)
    else:
        raise ValueError(coloring_property+" for coloring type not defined")
    if data_processed is not None:
        if coloring_property == "Missingness":
            column_values_processed = np.sum(~np.isnan(data_processed.values),axis=0)
        elif coloring_property=="Mean_Intensities":
            column_values_processed = np.nanmean(data_processed.values,axis=0)
        elif coloring_property=="Median_Intensities":
            column_values_processed = np.nanmedian(data_processed.values,axis=0)
        elif coloring_property=="Sum_Intensities":
            column_values_processed = np.log(np.nansum(base**(data_processed.values),axis=0))/np.log(base)
        return column_values_origin, column_values_processed
    return column_values_origin

def get_property(data, indices, metadata=None, column_metadata=None, plot_by="Correlation", data_type="continuous", shuffle=False, seed=42, take_abs = False, min_samples=2):
    """
    Calculate property for each sample.

    Parameters
    ----------
    data : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples).
    indices : list
        Indices of samples to consider for calculation.
    metadata : pd.DataFrame, optional
	    Metadata for each sample with the same indices as the columns of data_origin. The default is None.
    column_metadata : str, optional
        Column in Metadata containing the covariate for correlation (plot_by == "Correlation" or plot_by == "Corr_pvalue"). The default is None.
    plot_by : str, optional
        Property for which data should be calculated. Options: "Missingness", "Mean", "Median", "CV", "norm_MAE", "norm_MSE", "Correlation", "Corr_pvalue". The default is "Correlation".
    data_type : str, optional
        If plot_by == "Correlation" or plot_by == "Corr_pvalue", type of data that is correlated. Options: "continuous", "ordinal", "nominal". The default is "continuous".
    shuffle : bool, optional
        If True, the metadata are shuffled before correlation. The default is False.
    seed : int, optional
        Seed for shuffeling the data. The default is 42.
    take_abs : bool, optional
        If True, take the absolute value for correlation. The default is False.
    min_samples : int, optional
        Minimal number of samples with defined value. The default is 2.

    Returns
    -------
    property_peptide : pd. Series
        Calculated property per sample.
    """
    # Avoid errors if subset chosen is too small -> test
    #data_index = [index for index in data.index if sum(~data.loc[index,indices].isna()) >= min_samples]
    if min_samples < 2:
        min_samples = 2
    data_index = data.loc[((~data.loc[:,indices].isna()).T.sum() >= min_samples),:].index

    if plot_by == "Missingness":
        property_peptide = data.loc[data_index,indices].isna().T.sum()
    elif plot_by == "Mean":
        property_peptide = data.loc[data_index,indices].T.mean()
    elif plot_by == "Median":
        property_peptide = data.loc[data_index,indices].T.median()
    elif plot_by == "CV":
        property_peptide = pd.Series(index=data_index, data=[get_CV_log(data.loc[index,indices]) for index in data_index])
    elif plot_by == "norm_MAE":
        property_peptide = pd.Series(index=data_index, data = [mean_absolute_error(data.loc[index,indices].fillna(data.loc[index,indices].mean()),data.loc[index,indices].mean()*np.ones(len(indices)))/data.loc[index,indices].mean() for index in data_index])
    elif plot_by == "norm_MSE":
        property_peptide = pd.Series(index=data_index, data = [mean_squared_error(data.loc[index,indices].fillna(data.loc[index,indices].mean()),data.loc[index,indices].mean()*np.ones(len(indices)))/data.loc[index,indices].mean() for index in data_index])
    elif plot_by in ["Correlation","Corr_pvalue"]:
        if metadata is None or not column_metadata:
            raise ValueError("Metadata or Metadata column is not defined.")
        masks = ~np.isnan(data.loc[data_index,indices].values)
        meta_values = metadata.loc[indices, column_metadata].to_numpy()
        metadata_values_ = [meta_values[mask] for mask in masks]
        #metadata_values_ = [metadata.loc[data.loc[index,indices].dropna().index,column_metadata].values for index in data_index]
        if shuffle:
            np.random.seed(seed)
            for i, item in enumerate(metadata_values_):
                shuffled = np.array(item, copy=True)
                np.random.shuffle(shuffled)
                metadata_values_[i] = shuffled
                
        stats = []
        for ind, row in enumerate( data.loc[data_index,indices].values):
            mask = masks[ind]
        
            if data_type == "continuous":
                stats.append(pearsonr(metadata_values_[ind],row[mask]))
                #stats = [pearsonr(metadata_values_[ind],data.loc[index,data.loc[index,indices].dropna().index]) for ind,index in enumerate(data_index)]
            elif data_type == "ordinal":
                stats.append(spearmanr(metadata_values_[ind],row[mask]))
                #stats = [spearmanr(metadata_values_[ind],data.loc[index,data.loc[index,indices].dropna().index]) for ind,index in enumerate(data_index)]
            elif data_type == "nominal":
                stats.append(pointbiserialr(metadata_values_[ind],row[mask]))
                #stats = [pointbiserialr(metadata_values_[ind],data.loc[index,data.loc[index,indices].dropna().index]) for ind,index in enumerate(data_index)]
            else:
                raise ValueError("data_type not defined.")

        
        if plot_by == "Correlation":
            if take_abs:
                property_peptide = pd.Series(index=data_index, data=[abs(el[0]) for el in stats])
            else:
                property_peptide = pd.Series(index=data_index, data=[el[0] for el in stats])
        elif plot_by == "Corr_pvalue":
            property_peptide = pd.Series(index=data_index, data=[-1*np.log10(el[1]) for el in stats])
    return property_peptide

def set_aspect_equal(ax):
    """
    Set the aspect of the axis equal.

    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes
        Axes for which the action should be performed.
    """
    
    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()
    aspect = abs((x_lim[1]-x_lim[0])/(y_lim[1]-y_lim[0]))
    ax.set_aspect(aspect)
    return

def get_confidence_ellipse(x, y, confidence=0.95, **kwargs):
    """
    Calculation of confidence ellipses.

    Parameters
    ----------
    x : np.ndarray
        x-coordinate for the ellipse estimation.
    y : np.ndarray
        y-coordinate for the ellipse estimation.
    confidence : str, optional
        Confidence intervall used. The default is 0.95.

    Returns
    -------
    angle : float
        Angle of the ellipse.
    width : float
        Width of the ellipse.
    height : float
        Height of the ellipse.
    """

    # Covariance matrix
    cov = np.cov(x, y)

    # Eigen decomposition
    eigvals, eigvecs = np.linalg.eigh(cov)

    # Sort eigenvalues/vectors
    order = eigvals.argsort()[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    # Chi-square quantile
    chi2_val = chi2.ppf(confidence, df=2)
    #print(chi2_val)

    # Ellipse radii
    width = 2 * np.sqrt(eigvals[0] * chi2_val)
    height = 2 * np.sqrt(eigvals[1] * chi2_val)

    # Rotation angle in degrees
    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
    
    # Get axes extends
    #max_x = np.mean(x)+np.sqrt(chi2_val*cov[0,0])
    #min_x = np.mean(x)-np.sqrt(chi2_val*cov[0,0])

    #max_y = np.mean(y)+np.sqrt(chi2_val*cov[1,1])
    #min_y = np.mean(y)-np.sqrt(chi2_val*cov[1,1])

    return angle, width, height, #min_x, max_x, min_y, max_y

def get_vectors(metadata, data_transformed, data_columns, n_components = 4, min_samples = 3, max_entries = 10, max_entries_forced = 30, categorical_properties = None, verbose = True):
    """
    Calculation of the correlation vectors of a single transformed dataset.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties).
    data_transformed : np.ndarray
        Sample coordinates of the principal components (Indices: Samples, Columns: Components).
    data_columns : list
        Column names of the data, in the row order of data_transformed.
    n_components : int, optional
        First (highest Eigenvalue) n_components the properties are projected on. The default is 4.
    min_samples : int, optional
        Minimum number of samples with a finite value in both the property and the component, below which the coordinate is set to NaN. The default is 3.
    max_entries : int, optional
        Maximum number of levels of a categorical property, above which the property is dropped. The default is 10.
    max_entries_forced : int, optional
        Maximum number of levels of a property matched by categorical_properties. The default is 30.
    categorical_properties : tuple of str, optional
        Substrings marking properties that are one-hot encoded regardless of their dtype. The default is None.
    verbose : bool, optional
        If True, the dropped and the re-encoded properties are printed. The default is False.

    Returns
    -------
    vectors : pd.DataFrame
        Coordinates of each property (Indices: Properties/Levels, Columns: PC_1 ... PC_n) in the range [-1,1].

    """
    metadata_, data_transformed, components = get_aligned_metadata(metadata, data_transformed, data_columns, n_components)
    design = get_design_matrix(metadata_, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, verbose = verbose)

    vectors = pd.DataFrame(np.nan, index = design.columns, columns = components, dtype = float)
    for property_name in design.columns:
        values = design[property_name].to_numpy(dtype = float)
        for i, component in enumerate(components):
            coordinates = data_transformed[:,i]
            mask = np.isfinite(values) & np.isfinite(coordinates)
            if mask.sum() < min_samples:
                continue
            if (values[mask].std() < 1e-12) or (coordinates[mask].std() < 1e-12):
                continue
            vectors.loc[property_name,component] = float(np.corrcoef(values[mask], coordinates[mask])[0,1])
    return vectors

def get_eta2(metadata, data_transformed, data_columns, n_components = 4, min_samples = 3, max_entries = 10, max_entries_forced = 30, categorical_properties = None, adjusted = False, verbose = True):
    """
    Calculation of the squared correlation ratio of a single transformed dataset.

    All levels of a property are fitted jointly, so that the returned value is the fraction of the variance of a component explained by that property.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties).
    data_transformed : np.ndarray
        Sample coordinates of the principal components (Indices: Samples, Columns: Components).
    data_columns : list
        Column names of the data, in the row order of data_transformed.
    n_components : int, optional
        First (highest Eigenvalue) n_components the properties are projected on. The default is 4.
    min_samples : int, optional
        Minimum number of samples with a finite value in both the property and the component, below which the value is set to NaN. The default is 3.
    max_entries : int, optional
        Maximum number of levels of a categorical property, above which the property is dropped. The default is 10.
    max_entries_forced : int, optional
        Maximum number of levels of a property matched by categorical_properties. The default is 30.
    categorical_properties : tuple of str, optional
        Substrings marking properties that are one-hot encoded regardless of their dtype. The default is None.
    adjusted : bool, optional
        If True, eta2 is corrected for the number of levels of a property. The default is False.
    verbose : bool, optional
        If True, the dropped and the re-encoded properties are printed. The default is False.

    Returns
    -------
    eta2 : pd.DataFrame
        Squared correlation ratio of each property (Indices: Properties, Columns: PC_1 ... PC_n) in the range [0,1].

    """
    metadata_, data_transformed, components = get_aligned_metadata(metadata, data_transformed, data_columns, n_components)
    design = get_design_matrix(metadata_, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, verbose = verbose)

    blocks = {}
    for property_name in metadata_.columns:
        levels = [entry for entry in design.columns if (entry == property_name) or str(entry).startswith("%s=" % property_name)]
        if len(levels) > 0:
            blocks[property_name] = levels

    eta2 = pd.DataFrame(np.nan, index = list(blocks.keys()), columns = components, dtype = float)
    for property_name, levels in blocks.items():
        block = design[levels].to_numpy(dtype = float)
        for i, component in enumerate(components):
            coordinates = data_transformed[:,i]
            mask = np.isfinite(coordinates) & np.all(np.isfinite(block), axis = 1)
            if mask.sum() < min_samples:
                continue
            target = coordinates[mask]
            total = float(((target - target.mean())**2).sum())
            if total < 1e-12:
                continue
            predictors = np.column_stack([np.ones(int(mask.sum())), block[mask,:]])
            coefficients = np.linalg.lstsq(predictors, target, rcond = None)[0]
            residual = float(((target - predictors.dot(coefficients))**2).sum())
            value = 1.0 - (residual/total)
            degrees = int(np.linalg.matrix_rank(predictors)) - 1
            if adjusted and ((mask.sum() - degrees - 1) > 0):
                value = 1.0 - ((1.0 - value)*(mask.sum() - 1)/(mask.sum() - degrees - 1))
            elif not adjusted:
                value = min(max(value, 0.0), 1.0)
            eta2.loc[property_name,component] = value
    return eta2

def get_aligned_metadata(metadata, data_transformed, data_columns, n_components):
    """
    Alignment of the metadata with the samples of a transformed dataset.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties).
    data_transformed : np.ndarray
        Sample coordinates of the principal components (Indices: Samples, Columns: Components).
    data_columns : list
        Column names of the data, in the row order of data_transformed.
    n_components : int
        Requested number of principal components.

    Returns
    -------
    metadata_ : pd.DataFrame
        Metadata restricted to the samples of the dataset and ordered like data_transformed.
    data_transformed : np.ndarray
        Sample coordinates restricted to the first n_components components.
    components : list of str
        Names of the components ("PC_1" ... "PC_n").

    """
    data_columns = list(data_columns)
    data_transformed = np.asarray(data_transformed, dtype = float)
    if data_transformed.ndim != 2:
        raise ValueError("data_transformed has to be two-dimensional (Indices: Samples, data_columns: Components).")
    if data_transformed.shape[0] != len(data_columns):
        raise ValueError("data_transformed (%i samples) and data_columns (%i samples) do not match." % (data_transformed.shape[0], len(data_columns)))
    missing = [entry for entry in data_columns if entry not in metadata.index]
    if len(missing) > 0:
        raise ValueError("%i samples are not contained in the metadata, e.g. %s." % (len(missing), missing[0]))

    metadata_ = metadata.loc[data_columns]
    n_components = min(n_components, data_transformed.shape[1])
    components = ["PC_%i" % (i+1) for i in range(n_components)]
    return metadata_, data_transformed[:,:n_components], components

def get_design_matrix(metadata, max_entries = 10, max_entries_forced = 30, categorical_properties = None, verbose = False):
    """
    Transformation of the metadata into a purely numeric design matrix.

    Numeric and boolean properties are kept as they are, every other property is one-hot encoded into one indicator column per level ("Property=Level"). 
    Missing values are kept as NaN, also in the indicator columns, so that samples without a label are excluded instead of being counted as a zero. 
    Properties with exactly two levels are encoded by a single indicator column, since the second one would only mirror the first.

    Numeric properties whose name matches categorical_properties are one-hot encoded as well. 
    This covers label-encoded categories such as a batch numbered 1,2,3, whose numbering carries no order and whose effect would be invisible in a linear orrelation unless it happens to be monotone in the numbering.

    Categorical properties with a single level carry no information and categorical properties with too many levels (identifiers, free text) are too fragmented to be interpretable, both are dropped. T
    he upper limit is max_entries, or the usually more generous max_entries_forced for the properties matched by categorical_properties, whose number of levels is not known beforehand.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties).
    max_entries : int, optional
        Maximum number of levels (excluding missing values) of a categorical property, above which the property is dropped. The default is 10.
    max_entries_forced : int, optional
        Maximum number of levels of a property matched by categorical_properties, above which the property is dropped. Set generously, since the number of batches of a study is not known beforehand, but low enough to still catch free text. The default is 30.
    categorical_properties : tuple of str, optional
        Substrings, matched case-insensitively against the column names, marking properties that are one-hot encoded regardless of their dtype and limited by max_entries_forced instead of max_entries. An exact column name works as well, since it is a substring of itself. The default is None.
    verbose : bool, optional
        If True, the dropped properties and the numeric properties encoded as categories are printed. The default is False.

    Returns
    -------
    design : pd.DataFrame
        Numeric representation of the metadata (Indices: Samples, Columns: Properties/Levels).

    """
    categorical_properties = [] if categorical_properties is None else [str(entry).lower() for entry in categorical_properties]

    blocks = []
    dropped = {}
    encoded = []
    for property_name in metadata.columns:
        values = metadata[property_name]
        is_forced = any([entry in str(property_name).lower() for entry in categorical_properties])
        is_numeric = pd.api.types.is_bool_dtype(values) or pd.api.types.is_numeric_dtype(values)

        if is_numeric and not is_forced:
            blocks.append(values.astype(float).rename(property_name))
            continue
        if is_numeric:
            encoded.append(property_name)
            if np.all(np.mod(values.dropna().to_numpy(dtype = float), 1) == 0):
                values = values.astype("Int64")

        n_entries = values.dropna().nunique()
        limit = max_entries_forced if is_forced else max_entries
        if (n_entries < 2) or (n_entries > limit) or (n_entries >= len(metadata)):
            dropped[property_name] = n_entries
            continue

        dummies = pd.get_dummies(values, prefix = property_name, prefix_sep = "=", dtype = float)
        if n_entries == 2:
            dummies = dummies.iloc[:,[1]]
        dummies.loc[values.isna(),:] = np.nan
        blocks.append(dummies)

    if verbose and (len(encoded) > 0):
        print("Encoded %i numeric properties as categories: %s" % (len(encoded), ", ".join([str(entry) for entry in encoded])))
    if verbose and (len(dropped) > 0):
        print("Dropped %i categorical properties (allowed: 2-%i levels, 2-%i for %s): %s" % (len(dropped), max_entries, max_entries_forced, "/".join(categorical_properties) if len(categorical_properties) > 0 else "none", ", ".join(["%s (%i levels)" % (key,value) for key,value in dropped.items()])))
    if len(blocks) == 0:
        raise ValueError("metadata does not contain any property that can be projected.")
    return pd.concat(blocks, axis = 1)

def get_metadata_property_colors(entries, metadata, cmap = "tab10"):
    """
    Assignment of one color per metadata property to the rows of a vector or eta2 table.

    Parameters
    ----------
    entries : list
        Row names of the table ("Property" or "Property=Level").
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties).
    cmap : str, optional
        Colormap the colors are taken from. The default is "tab10".

    Returns
    -------
    properties : list
        Property each entry belongs to, in the order of entries.
    labels : list
        Level of each entry, or the property itself for a property without levels.
    colors : dict
        Color for each property.

    """
    properties = []
    labels = []
    for entry in entries:
        name = str(entry)
        for property_name in metadata.columns:
            if (entry == property_name) or name.startswith("%s=" % property_name):
                name = str(property_name)
                break
        properties.append(name)
        labels.append(str(entry)[len(name)+1:] if str(entry) != name else name)

    palette = plt.get_cmap(cmap)
    unique = list(dict.fromkeys(properties))
    colors = {name: palette(index % palette.N) for index, name in enumerate(unique)}
    return properties, labels, colors

def get_visible_entries(table, i, j, top_n = None, noise_level = None):
    """
    Selection of the entries of a vector or eta2 table drawn within the plane of two principal components.

    The entries are ranked by their distance from the origin within that plane, which is the length of the vector for get_metadata_vectors and the distance within the relationship square for get_metadata_eta2. 
    In both cases the squared distance is the fraction of the variance of the property that the two components account for, so a property that is barely represented within the plotted plane is ranked last.
    Entries without a finite value in one of the two components are never drawn.

    Parameters
    ----------
    table : pd.DataFrame
        Coordinates of each property (Indices: Properties/Levels, Columns: PC_1 ... PC_n).
    i : int
        Component plotted along the x-axis.
    j : int
        Component plotted along the y-axis.
    top_n : int, optional
        If set, only the top_n entries with the highest distance are returned. The default is None.
    noise_level : float, optional
        If set, entries with a distance below noise_level are not returned. The default is None.

    Returns
    -------
    entries : list
        Rows to be drawn, ordered by decreasing distance from the origin.

    """
    x = table["PC_%i" % (i+1)]
    y = table["PC_%i" % (j+1)]
    distances = np.sqrt(x**2 + y**2)
    distances = distances[np.isfinite(x) & np.isfinite(y)]
    if noise_level is not None:
        distances = distances[distances >= noise_level]
    distances = distances.sort_values(ascending = False)
    if top_n is not None:
        distances = distances.iloc[:top_n]
    return list(distances.index)

def get_axis_limits(table, n_components, axis_offset = 0.1):
    """
    Calculation of one axis limit per principal component, shared by the whole row and column of a dataset.

    Parameters
    ----------
    table : pd.DataFrame
        Coordinates of each property (Indices: Properties, Columns: PC_1 ... PC_n).
    n_components : int
        Number of components shown in the grid.
    axis_offset : float, optional
        Fraction of the span added as a margin, so that a property lying on the origin is not cut by the axis. The default is 0.1.

    Returns
    -------
    limits : list of tuple of float
        Lower and upper limit of each component, in the order PC_1 ... PC_n.

    """
    components = ["PC_%i" % (index+1) for index in range(n_components)]
    values = table[components].to_numpy(dtype = float)
    values = np.where(np.isfinite(values), values, np.nan)
    if not np.isfinite(values).any():
        return [(-0.05,1.05) for component in components]

    lower = min(0.0, float(np.nanmin(values)))
    span = float(np.nanmax(values)) - lower
    if span < 1e-12:
        span = 1.0
    margin = span*axis_offset/2.0
    return [(lower-margin, lower+span+margin) for component in components]

def merge_duplicate_entries(table, separator = " / ", decimals = 12):
    """
    Merging of properties that have identical coordinates into a single entry.
 
    Parameters
    ----------
    table : pd.DataFrame
        Coordinates of each property (Indices: Properties, Columns: PC_1 ... PC_n).
    separator : str, optional
        String the merged property names are joined with. The default is " / ".
    decimals : int, optional
        Number of decimals the coordinates are rounded to before they are compared. The default is 12.
 
    Returns
    -------
    merged : pd.DataFrame
        Table with one row per unique set of coordinates, in the order of the input.
 
    """
    keys = {}
    for property_name in table.index:
        key = tuple(table.loc[property_name].round(decimals).fillna(np.inf))
        keys.setdefault(key, []).append(str(property_name))
 
    merged = table.loc[[names[0] for names in keys.values()]].copy()
    merged.index = [separator.join(names) for names in keys.values()]
    return merged

### For all modules to check how either plot changes, adjust labels for summarization
def plot_qa_box_int_and_number(data_origin, metadata, data_processed=None, sample_type=None, sample_type_column=None, size = (10, 10), fontsize = 15, ylabel_int = "# Precursors", y_label_box = "log2(Precursor quantity)", x_label = "Samples", step = "", outfile = "Intensity_distribution_boxplot.png", column_color = None, use_metadata_coloring = False, cmap=None, continuous=False):
    """
    Plotting of top: Number of measured IDs per sample and bottom: Boxplot of the intensities (25, 50, 75%) with whiskers extending max 1.5*IQR (IQR: 25-75%).

    Parameters
    ----------
    data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame
	    Metadata for each sample with the same indices as the columns of data_origin.
	data_processed : pd.DataFrame, optional
	    Data (Indices: Precursors, Columns: Samples) after the preprocessing step.
    sample_type : str, optional
	    Sample type in Metadata[sample_type_column] for which indices should be extracted. The default is None.
    sample_type_column : str, optional
	    Column in Metadata that distinguishes sample type. If None, samples will not be distinguished. The default is None.
    size : tuple of float, optional
        Size of the plot. The default is (10, 10).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    ylabel_int : str, optional
        y-axis label (top). The default is "# Precursors".
    y_label_box : str, optional
        y-axis label (bottom). The default is "log2(Precursor quantity)".
    x_label : str, optional
        x-axis label. The default is "Samples".
    step : str, optional
        Preprocessing step (used for the title). The default is "".
    outfile : str, optional
        File path for the figure. The default is "Intensity_distribution_boxplot.png".
    column_color : str, optional
        If use_metadata_coloring == True, Metadata column used for coloring. The default is None.
    use_metadata_coloring : bool, optional
        If True, Metadata are used for coloring. The default is False.
    cmap : str, optional
        Colormap used for coloring. The default is "viridis".
    continuous : bool, optional
        If True continuous data are assumed for coloring e.g. Order. If False discrete data are assumed e.g. Batch. The default is False.
    """
    
    # Setup
    if not cmap:
        if continuous:
            cmap = "viridis"
        else:
            cmap = "Spectral"
    
    fig = plt.figure(figsize=size)
    if data_processed is not None:
        gs = GridSpec(2,2,height_ratios = [0.25,0.75], width_ratios = [0.5,0.5], hspace = 0.05, wspace = 0.05)
    else:
        gs = GridSpec(2,1,height_ratios = [0.25,0.75], hspace= 0.05)

    ax_numbers_origin = fig.add_subplot(gs[0,0])
    ax_boxplot_origin = fig.add_subplot(gs[1,0])
    
    if data_processed is not None:
        ax_numbers_processed = fig.add_subplot(gs[0,1])
        ax_boxplot_processed = fig.add_subplot(gs[1,1])
    
    # sample_type selection
    if data_processed is not None:
        indices_origin, indices_processed = get_indices_by_type(data_origin, metadata, data_processed, sample_type_column, sample_type)
    else:
        indices_origin = get_indices_by_type(data_origin, metadata, data_processed, sample_type_column, sample_type)
    
    # metadata coloring
    if use_metadata_coloring:
        if data_processed is not None:
            colors_origin, colors_processed, norm, values_dict = get_colors_by_metadata(metadata, column_color, indices_origin, indices_processed=indices_processed, cmap=cmap, continuous=continuous)
        else:
            colors_origin, norm, values_dict = get_colors_by_metadata(metadata, column_color, indices_origin, indices_processed=None, cmap=cmap, continuous=continuous)
    
    # Plotting
    y_origin = [sum(~np.isnan(data_origin[column].values)) for column in indices_origin]
    if use_metadata_coloring:
        ax_numbers_origin.scatter(np.arange(len(indices_origin)),y_origin, c = colors_origin)
    else:
        ax_numbers_origin.scatter(np.arange(len(indices_origin)),y_origin)
    ax_numbers_origin.set_xticks(np.arange(len(indices_origin)))
    ax_numbers_origin.set_xticklabels([])
    ax_numbers_origin.set_ylabel(ylabel_int, fontsize=fontsize)
    ax_numbers_origin.set_axisbelow(True)
    ax_numbers_origin.grid(axis="y")
    ax_numbers_origin.tick_params(labelsize=fontsize)
    
    if data_processed is not None:
        y_processed = [sum(~np.isnan(data_processed[column].values)) for column in indices_processed]
        if use_metadata_coloring:
            ax_numbers_processed.scatter(np.arange(len(indices_processed)),y_processed, c=colors_processed)
        else:
            ax_numbers_processed.scatter(np.arange(len(indices_processed)),y_processed)
        ax_numbers_processed.set_xticks(np.arange(len(indices_processed)))
        ax_numbers_processed.set_xticklabels([])
        ax_numbers_processed.set_yticklabels([])
        ax_numbers_processed.set_axisbelow(True)
        ax_numbers_processed.grid(axis="y")
        ax_numbers_processed.tick_params(labelsize=fontsize)
        
    if data_processed is not None:
        ax_numbers_origin.set_title("Before "+step, fontsize=fontsize)
        ax_numbers_processed.set_title("After "+step, fontsize=fontsize)
    
    if data_processed is not None:
        y_lim = (np.nanmin([np.nanmin(y_origin),np.nanmin(y_processed)])*0.95,np.nanmax([np.nanmax(y_origin),np.nanmax(y_processed)])*1.05)
    else:
        y_lim = (np.nanmin(y_origin)*0.95,np.nanmax(y_origin)*1.05)
        
    ax_numbers_origin.set_ylim(y_lim)
    ax_numbers_origin.set_xlim(-1,len(indices_origin))
    if data_processed is not None:
        ax_numbers_processed.set_ylim(y_lim)
        ax_numbers_processed.set_xlim(-1,len(indices_processed))

    data_origin_filtered = [data_origin[column].copy().dropna().values for column in indices_origin]
    bplot = ax_boxplot_origin.boxplot(data_origin_filtered, patch_artist=True, medianprops=dict(color="w"))
    
    if use_metadata_coloring:
        for patch, color in zip(bplot['boxes'], colors_origin):
            patch.set_facecolor(color)
    
    ax_boxplot_origin.set_ylabel(y_label_box, fontsize=fontsize)
    ax_boxplot_origin.set_xticklabels([])
    ax_boxplot_origin.set_xlabel(x_label, fontsize=fontsize)
    ax_boxplot_origin.tick_params(labelsize=fontsize)
    ax_boxplot_origin.set_axisbelow(True)
    ax_boxplot_origin.grid(axis="y")

    if data_processed is not None:
        data_processed_filtered = [data_processed[column].copy().dropna().values for column in indices_processed]
        bplot = ax_boxplot_processed.boxplot(data_processed_filtered, patch_artist=True, medianprops=dict(color="w"))
        
        if use_metadata_coloring:
            for patch, color in zip(bplot['boxes'], colors_processed):
                patch.set_facecolor(color)
        
        ax_boxplot_processed.set_xticklabels([])
        ax_boxplot_processed.set_xlabel(x_label, fontsize=fontsize)
        ax_boxplot_processed.set_yticklabels([])
        ax_boxplot_processed.tick_params(labelsize=fontsize)
        ax_boxplot_processed.set_axisbelow(True)
        ax_boxplot_processed.grid(axis="y")
        
    if data_processed is not None:
        y_lim = (np.min([data_processed[indices_processed].replace(-1*np.inf,np.nan).min().min(),data_origin[indices_origin].replace(-1*np.inf,np.nan).min().min()])*0.95,
                 np.max([data_processed[indices_processed].replace(+1*np.inf,np.nan).max().max(),data_origin[indices_origin].replace(+1*np.inf,np.nan).max().max()])*1.05)
    else:
        y_lim = (data_origin[indices_origin].min().min()*0.95,data_origin[indices_origin].max().max()*1.05)

    ax_boxplot_origin.set_ylim(y_lim)
    ax_boxplot_origin.set_xlim(0,len(indices_origin)+1)
    if data_processed is not None:
        ax_boxplot_processed.set_ylim(y_lim)
        ax_boxplot_processed.set_xlim(0,len(indices_processed)+1)

    plt.savefig(outfile, bbox_inches="tight")
    return

### For Filtering_long and Filtering_wide
def plot_filter_mask(mask, mask_description, mask_description_base=None, use_mask_values=False, size=(5,5), fontsize = 15, fontsize_filters=10, outfile="Mask_description_barplot.png"):
    """
    Plotting of the filter mask

    Parameters
    ----------
    mask : pd.DataFrame
        Filter mask.
    mask_description : pd.DataFrame
        Unique identifier for each filter.
    mask_description_base : pd.DataFrame, optional
        Unique identifiers for each filter of previous steps. The default is None.
    use_mask_values : bool, optional
        If True, mask values extracted from the mask instead of the mask description and used for plotting. The default is False.
    size : tuple of float, optional
        Size of the plot. The default is (5,5).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    fontsize_filters : int, optional
        Font size of the x-tick labels. The default is 10.
    outfile : str, optional
        File path for the figure. The default is "Mask_description_barplot.png".
    """    
    
    if use_mask_values:
        mask.fillna(0,inplace=True)
        mask_values = list(set(pd.melt(mask)["value"].values)^set([0]))
        mask_values.sort()
    else:
        description_dict = {mask_description.loc[index,"value"]:index for index in mask_description.index}
        if mask_description_base is not None:
            mask_values = [value for value in list(description_dict.keys()) if value not in [0,1,99] and value not in list(mask_description_base["value"].values)]
        else:
            mask_values = [value for value in list(description_dict.keys()) if value not in [0,1,99]]
    
    heigths = [(mask==value).sum().sum() for value in mask_values]

    fig,ax = plt.subplots()
    fig.set_size_inches(size)
    
    ax.bar(np.arange(len(mask_values)),heigths)
    
    ax.set_xticks(np.arange(len(mask_values)),heigths)
    if use_mask_values:
        ax.set_xticklabels(mask_values, fontsize = fontsize_filters, rotation=45, ha="right")
    else:
        ax.set_xticklabels(np.vectorize(description_dict.get)(mask_values),fontsize = fontsize_filters, rotation=45, ha="right")
    ax.tick_params(axis="y", labelsize = fontsize)
    ax.set_ylabel("Occurence", fontsize = fontsize)
    ax.set_title("Dropouts", fontsize = fontsize)
    ax.set_yscale("log")
    ax.set_axisbelow(True)
    ax.grid(axis="y",which="major")
        
    plt.savefig(outfile, bbox_inches="tight")
    return

# For Imputation, Normalization, Correction, Summarization
def plot_CV_violin(data_origin, metadata, data_processed = None, group_columns = [], size=(1,5), fontsize = 15, fontsize_title = 15, calc_widths=False, step="", outfile="CV_violinplot.png", base=2, scale_max=0.95, y_lim=None, write_to_disk=False):
    """
    Plotting of the coefficient of variation (CV) as a violin plot for different sample types.

    Parameters
    ----------
	data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame
	    Metadata for each sample with the same indices as the columns of data_origin.
	data_processed : pd.DataFrame, optional
	    Data (Indices: Precursors, Columns: Samples) after the preprocessing step.
    group_columns : list of str, optional
        Columns in Matedata that distinguish different sample types. This information can also contain other information e.g. "Batch". The default is [].
    size : tuple of float, optional
        Size of the plot. The default is (1,5).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    fontsize_title : int, optional
        Fontsize of the title. The default is 15.
    calc_widths : bool, optional
        If True scales violins by maximal width. The default is True.
    step : str, optional
        Preprocessing step (used for the title). The default is "".
    outfile : str, optional
        File path for the figure. The default is "CV_violinplot.png".
    base : float, optional
        Base of the logarithmic transformation. The default is 2.
    scale_max : float, optional
        Maximal width for violin plot. The default is 0.95.
    y_lim : tuple of float, optional
        y-limit for which the CVs should be evaluated (to mask outliers for better visuals ONLY). If None, no filter is applied. The default is None.
    write_to_disk : bool, optional
        If true saves CV to disk as outfile[:-4]+"_CVs_origin.tsv". The default is False.
        
    Returns
    -------
    cross_border_origin : bool
        Only returned if y_lim is not None. If True is returned at least one CV crossed max(y_lim) in original data.
    cross_border_processed : bool
        Only returned if y_lim is not None. If True is returned at least one CV crossed max(y_lim) in processed data.
    
    """
    
    if data_processed is not None:
        group_dict_comp = {str(group_column)+"__"+str(group):[index for index in metadata[metadata[group_column]==group].index if index in data_origin.columns and index in data_processed.columns] for group_column in group_columns for group in np.unique(metadata[group_column].values)}
    else:
        group_dict_comp = {str(group_column)+"__"+str(group):[index for index in metadata[metadata[group_column]==group].index if index in data_origin.columns] for group_column in group_columns for group in np.unique(metadata[group_column].values)}
    group_dict = {}
    
    keys = list(group_dict_comp.keys())
    for ind,key in enumerate(keys):
        for key_2 in keys[ind+1:]:
            if key.split("__")[0] != key_2.split("__")[0]:
                group_dict.update({key+"___"+key_2:list(set(group_dict_comp[key])&set(group_dict_comp[key_2]))})
    if not group_dict:
        group_dict = group_dict_comp
    if not group_dict:
        group_dict = {"All":list(data_origin.columns)}  
        
    #cvs_lognormal_origin = [[el for el in [get_CV_log(data_origin.loc[index,[column for column in group_dict[group] if column in data_origin.columns]], base=base) for index in data_origin.index]] for group in group_dict.keys()]
    cvs_lognormal_origin = [get_CV_log_array(data_origin.loc[:,[column for column in group_dict[group] if column in data_origin.columns]], base=base) for group in group_dict.keys()]

    if write_to_disk:
        pd.DataFrame(columns=data_origin.index,index=group_dict.keys(),data=cvs_lognormal_origin).T.to_csv(outfile[:-4]+"_CVs_origin.tsv",sep="\t")

    cvs_lognormal_origin = [[el for el in item if not np.isnan(el)] for item in cvs_lognormal_origin]
    cvs_lognormal_origin = [item if item or sum(item)>0 else [np.nan]*2 for item in cvs_lognormal_origin]
    
    if data_processed is not None:
        #cvs_lognormal_processed = [[el for el in [get_CV_log(data_processed.loc[index,[column for column in group_dict[group] if column in data_processed.columns]], base=base) for index in data_processed.index]] for group in group_dict.keys()]
        cvs_lognormal_processed = [get_CV_log_array(data_processed.loc[:,[column for column in group_dict[group] if column in data_processed.columns]], base=base) for group in group_dict.keys()]
        
        if write_to_disk:
            pd.DataFrame(columns=data_processed.index,index=group_dict.keys(),data=cvs_lognormal_processed).T.to_csv(outfile[:-4]+"_CVs_processed.tsv",sep="\t")
        
        cvs_lognormal_processed = [[el for el in item if not np.isnan(el)] for item in cvs_lognormal_processed]
        cvs_lognormal_processed = [item if item else [np.nan]*2 for item in cvs_lognormal_processed]
    
    if y_lim:
        cross_border_origin = any(True if el > y_lim[1] else False for el in [max(item) for item in cvs_lognormal_origin])
        cvs_lognormal_origin = [[el for el in item if el>=y_lim[0] and el<=y_lim[1]] for item in cvs_lognormal_origin]
        cvs_lognormal_origin = [[np.nan,np.nan] if not el else el for el in cvs_lognormal_origin]
        if data_processed is not None:
            cvs_lognormal_processed = [[el for el in item if el>=y_lim[0] and el<=y_lim[1]] for item in cvs_lognormal_processed]
            cvs_lognormal_processed = [[np.nan,np.nan] if not el else el for el in cvs_lognormal_processed]
            cross_border_processed = any(True if el > y_lim[1] else False for el in [max(item) for item in cvs_lognormal_processed])
        else:
            cross_border_processed = False
            
    if calc_widths:
        min_origin = np.nanmin([el for item in cvs_lognormal_origin for el in item])
        max_origin = np.nanmax([el for item in cvs_lognormal_origin for el in item])
        widths_origin = [np.nanmax(np.histogram(cvs_lognormal_origin[ind],bins=25, range=(min_origin,max_origin))[0]) for ind in range(len(cvs_lognormal_origin)) if cvs_lognormal_origin[ind]]
        if data_processed is not None:
            min_processed = np.nanmin([el for item in cvs_lognormal_processed for el in item])
            max_processed = np.nanmax([el for item in cvs_lognormal_processed for el in item])
            widths_processed = [np.nanmax(np.histogram(cvs_lognormal_processed[ind],bins=25, range=(min_processed,max_processed))[0]) for ind in range(len(cvs_lognormal_processed)) if cvs_lognormal_processed[ind]]
            
            scale = max(max(widths_origin),max(widths_processed))
            widths_origin = widths_origin/scale*scale_max
            widths_processed = widths_processed/scale*scale_max
        else:
            widths_origin = widths_origin/max(widths_origin)*scale_max
    else:
        widths_origin = scale_max
        widths_processed = scale_max
    
    size = ((size[0]*len(group_dict))+1,size[1])
    
    fig,ax = plt.subplots()
    fig.set_size_inches(size)
    
    bg_origin = ax.violinplot(cvs_lognormal_origin, side="low", showmedians=False, widths=widths_origin, showextrema=False)
    for pc in bg_origin['bodies']:
        pc.set_facecolor("w")
        pc.set_edgecolor("w")
        pc.set_alpha(1)
    
    core_origin = ax.violinplot(cvs_lognormal_origin, side="low", showextrema=True, showmedians=True, quantiles=[[0.25,0.75]]*len(cvs_lognormal_origin), widths=widths_origin)
    
    core_origin["cquantiles"].set_linestyle(":")
    core_origin["cquantiles"].set_linewidth(1)
    core_origin["cmedians"].set_linewidth(2.5)
    core_origin["cbars"].set_color("black")
    
    for pc in core_origin['bodies']:
        pc.set_facecolor('royalblue')
        pc.set_edgecolor('black')
    
    for part in ["cmaxes","cmins","cquantiles","cmedians"]:
        core_origin[part].set_color(["royalblue"]*len(cvs_lognormal_origin))
    
    if data_processed is not None:
        
        bg_processed = ax.violinplot(cvs_lognormal_processed, side="high", showmedians=False, widths=widths_processed, showextrema=False)
        for pc in bg_processed['bodies']:
            pc.set_facecolor("w")
            pc.set_edgecolor("w")
            pc.set_alpha(1)

        core_processed = ax.violinplot(cvs_lognormal_processed, side="high", showextrema=True, showmedians=True, quantiles=[[0.25,0.75]]*len(cvs_lognormal_processed), widths=widths_processed)
        
        core_processed["cquantiles"].set_linestyle(":")
        core_processed["cquantiles"].set_linewidth(1)
        core_processed["cmedians"].set_linewidth(2.5)
        core_processed["cbars"].set_color("black")
        
        for pc in core_processed['bodies']:
            pc.set_facecolor('orange')
            pc.set_edgecolor('black')
        
        for part in ["cmaxes","cmins","cquantiles","cmedians"]:
            core_processed[part].set_color(["orange"]*len(cvs_lognormal_processed))
    
        #ax.set_title("Before | After", fontsize = fontsize)
        fig.text(0.48, 0.9, "Before "+step, ha="right", va="bottom",color="royalblue", fontsize = fontsize_title)
        fig.text(0.5, 0.9, "|", ha="center", va="bottom", fontsize = fontsize_title)
        fig.text(0.52,0.9,"After "+step, ha="left", va="bottom",color="orange", fontsize = fontsize_title)
        
    ax.set_xticks(np.arange(1,len(group_dict)+1))
    ax.set_xticklabels([item.replace("___","\n").replace("__",":") for item in group_dict.keys()], fontsize = fontsize, rotation = 45, ha = "right")
    ax.set_ylabel("Coefficient of variation", fontsize = fontsize)
    ax.tick_params(axis="y", labelsize = fontsize)
    
    if y_lim:
        ax.set_ylim(y_lim)
    
    ax.set_axisbelow(True)
    ax.grid(axis="y")

    plt.savefig(outfile, bbox_inches="tight")
    if y_lim:
        return cross_border_origin, cross_border_processed
    else:
        return

# For Imputation, Normalization, Correction
### Coloring
# Corrections: columns used -> continuous scale
# Normalization: mean intensities and summed intensities (log) -> continuous scale
# Imputation: mean intensities and summed intensities (log) and missingness -> continuous scale
def plot_PCA(data_origin, metadata, data_processed = None, n_components = 4, size = (10,10), fontsize = 15, fontsize_title = 20, line_offset = 0.025, fontsize_label = 15, outfile = "PCA.png", continuous = True, column_color = None, use_metadata_coloring = False, coloring_property = "Missingness", base=2, cmap=None, plot_confidence_ellipse=False, confidence=0.95, set_axis_lim = False, axis_offset = 0.1, return_fitted_data = False):
    """
    Plotting of the first n_components principal components with different coloring options.

    Parameters
    ----------
    data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame
	    Metadata for each sample with the same indices as the columns of data_origin.
	data_processed : pd.DataFrame, optional
	    Data (Indices: Precursors, Columns: Samples) after the preprocessing step.
    n_components : int, optional
        First (highest Eigenvalue) n_components are plotted. The default is 4.
    size : tuple of float, optional
        Size of the plot. The default is (10,10).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    fontsize_title : int, optional
        Fontsize of the title. The default is 15.
    line_offset : float, optional
        Offset of diagonal line (needs to adjusted if a lot of components are plotted. The default is 0.025.
    fontsize_label : int, optional
        Font size of the colorbar label. The default is 15.
    outfile : str, optional
        File path for the figure. The default is "PCA.png".
    continuous : bool, optional
        If True continuous data are assumed for coloring e.g. Missingness. If False discrete data are assumed e.g. Batch. The default is True.
    column_color : str, optional
        If use_metadata_coloring == True, Metadata column used for coloring. The default is None.
    use_metadata_coloring : bool, optional
        If True, Metadata are used for coloring. The default is False.
    coloring_property : str, optional
        Property used for coloring. Options: "Missingness", Mean_Intensities", "Median_Intensities", Sum_Intensities" The default is "Missingness".
    base : float, optional
        Base of the logarithmic transformation. The default is 2.
    cmap : str, optional
        Colormap used for coloring. If None is selected, viridis is used as a default for continuous data and spectral for discrete data. The default is None.
    plot_confidence_ellipse : bool, optional
        If True, confidence ellipses are plotted. The default is False.
    confidence : float, optional
        Confidence interval for the ellipses. The default is 0.95.
    set_axis_lim : bool, optional
        Set axis-limits automatically (sometimes needed if plot_confidence_ellipse == True). The default is False.
    axis_offset : float, optional
        Axis offset for automatic axis-limits. The default is 0.1.
    return_fitted_data : bool, optional
        If True returns the fitted data for additional analysis
        
    Returns
    -------
    data_origin_transformed : np.ndarray, optional
        Fitted model for origin data
    data_processed_transformed : np.ndarray, optional
        Fitted model for processed data
    data_origin_columns : list, optional
        Column names of data_origin
    data_preocessed_columns : list, optional
        Column names of data_processed
    
    """
    
    ### Examples:
    ### Drift correction
    #continuous = True,
    #column_color = "Order",
    #use_metadata = False,
    
    ### Correction
    #continuous = False
    #column_color = "MS.Batch"
    #column_color = "Sample.Type"
        
    if not use_metadata_coloring:
        if not continuous:
            continuous = True
            print("No discrete coloring implemented without use_metadata_color==True. Used continuous=True instead.")
    
    if not cmap:
        if continuous:
            cmap = "viridis"
        else:
            cmap = "Spectral"
            
    # Coloring property
    data_origin_pca = data_origin.copy()
    if data_processed is not None:
        data_processed_pca = data_processed.copy()
    
    ### Implement check for metadata coloring
    
    data_origin_pca[data_origin_pca==0] = np.nan
    if data_processed is not None:
        data_processed_pca[data_processed_pca==0] = np.nan
        if not use_metadata_coloring:
            column_values_origin, column_values_processed = get_color_values_by_property(data_origin_pca, coloring_property, data_processed_pca, base)
    elif not use_metadata_coloring:
        column_values_origin = get_color_values_by_property(data_origin, coloring_property, base=base)
    
    data_origin_pca.dropna(inplace=True)
    data_origin_scaled = StandardScaler().fit_transform(data_origin_pca.T)
    if use_metadata_coloring:
        column_values_origin = metadata.loc[data_origin_pca.columns,column_color].values
        
    pca_origin = PCA()
    data_origin_transformed = pca_origin.fit_transform(data_origin_scaled)[:,:n_components]
    var_origin = pca_origin.explained_variance_ratio_[:n_components]
    
    if data_processed is not None:
        data_processed_pca.dropna(inplace=True)
        data_processed_scaled = StandardScaler().fit_transform(data_processed_pca.T)
        if use_metadata_coloring:
            column_values_processed = metadata.loc[data_processed_pca.columns,column_color].values
        
        pca_processed = PCA()
        data_processed_transformed = pca_processed.fit_transform(data_processed_scaled)[:,:n_components]
        var_processed = pca_processed.explained_variance_ratio_[:n_components]    
        
        colors_origin, colors_processed, norm, values_dict = get_colors_selection(column_values_origin, column_values_processed=column_values_processed, cmap=cmap, continuous=continuous)
    else:
        colors_origin, norm, values_dict = get_colors_selection(column_values_origin, column_values_processed=None, cmap=cmap, continuous=continuous)

    # Plotting
    fig = plt.figure(figsize=(n_components*2.5,n_components*2.5))
    gs = GridSpec(n_components,n_components, hspace= 0.1, wspace = 0.1)
    
    for i in range(n_components-1):
        for j in range(i+1,n_components):
            ax = fig.add_subplot(gs[j,i])  
            ax.scatter(data_origin_transformed[:,i],data_origin_transformed[:,j],alpha=0.5, c=colors_origin)  
        
            if not continuous and plot_confidence_ellipse:
                for color in set(colors_origin):
                    
                    group_indices = [index for index in range(len(data_origin_transformed)) if colors_origin[index]==color]
                    subset_i = data_origin_transformed[group_indices,i]
                    subset_j = data_origin_transformed[group_indices,j]
                    
                    angle, width, height = get_confidence_ellipse(subset_i,subset_j,confidence=confidence)
                    ax.add_patch(Ellipse(
                        (np.mean(subset_i), np.mean(subset_j)),
                        width=width,height=height,angle=angle,
                        edgecolor=color,facecolor="None",linewidth=2))
                    
            if set_axis_lim:
                x_min = np.min(data_origin_transformed[:,i])
                x_max = np.max(data_origin_transformed[:,i])
                x_offset = np.max((np.abs(x_min),x_max))*(axis_offset)
                
                y_min = np.min(data_origin_transformed[:,j])
                y_max = np.max(data_origin_transformed[:,j])
                y_offset = np.max((np.abs(y_min),y_max))*(axis_offset)
                
                ax.set_xlim(x_min-x_offset, x_max+x_offset)
                ax.set_ylim(y_min-y_offset, y_max+y_offset)
            
            ax.tick_params(labelsize = fontsize)
            set_aspect_equal(ax)
            if i != 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("PC "+str(j+1)+" ("+"{:.1f}".format(var_origin[j]*100)+"%)", fontsize = fontsize, labelpad = 5)
                ax.yaxis.get_major_locator().set_params(integer=True) 
            if j != n_components-1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("PC "+str(i+1)+" ("+"{:.1f}".format(var_origin[i]*100)+"%)", fontsize = fontsize, labelpad = 5)
                ax.xaxis.get_major_locator().set_params(integer=True) 
            
            ax.set_axisbelow(True)
            ax.grid()
            
            if data_processed is not None:
                ax = fig.add_subplot(gs[i,j])
                ax.scatter(data_processed_transformed[:,j],data_processed_transformed[:,i],alpha=0.5, c=colors_processed)
                
                if not continuous and plot_confidence_ellipse:
                    for color in set(colors_origin):
                        
                        group_indices = [index for index in range(len(data_processed_transformed)) if colors_origin[index]==color]
                        subset_i = data_processed_transformed[group_indices,i]
                        subset_j = data_processed_transformed[group_indices,j]
                        
                        angle, width, height = get_confidence_ellipse(subset_j,subset_i,confidence=confidence)
                        ax.add_patch(Ellipse(
                            (np.mean(subset_j), np.mean(subset_i)),
                            width=width,height=height,angle=angle,
                            edgecolor=color,facecolor="None",linewidth=2))
    
                if set_axis_lim:
                    x_min = np.min(data_processed_transformed[:,j])
                    x_max = np.max(data_processed_transformed[:,j])
                    x_offset = np.max((np.abs(x_min),x_max))*(axis_offset)

                    y_min = np.min(data_processed_transformed[:,i])
                    y_max = np.max(data_processed_transformed[:,i])
                    y_offset = np.max((np.abs(y_min),y_max))*(axis_offset)

                    ax.set_xlim(x_min-x_offset, x_max+x_offset)
                    ax.set_ylim(y_min-y_offset, y_max+y_offset)
                
                ax.tick_params(labelsize = fontsize, labeltop=True, labelright=True, labelbottom=False, labelleft = False)
                set_aspect_equal(ax)
                
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")
                ax.yaxis.get_major_locator().set_params(integer=True) 
                
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position("top")
                ax.xaxis.get_major_locator().set_params(integer=True) 
    
                if j != n_components-1:
                    ax.set_yticklabels([])
                else:
                    ax.set_ylabel("PC "+str(i+1)+" ("+"{:.1f}".format(var_processed[i]*100)+"%)", fontsize = fontsize, labelpad = 5)
                    
                if i != 0:
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel("PC "+str(j+1)+" ("+"{:.1f}".format(var_processed[j]*100)+"%)", fontsize = fontsize, labelpad = 5)
        
            ax.set_axisbelow(True)
            ax.grid()
    
    if data_processed is not None:
        s = 2.5*(np.log2(n_components)+1)/n_components
        s_min = 0.01*n_components
        s_max = 1 - s_min
        fig.text(1.5/(s*n_components),(s*n_components-1.5)/(s*n_components), "Before", rotation=315, va="center", ha="right",color="royalblue", fontsize = fontsize_title)
        fig.text(1.5/(s*n_components),(s*n_components-1.5)/(s*n_components), "After", rotation=315, va="center", ha="left",color="orange", fontsize = fontsize_title)
        line = plt.Line2D((s_min,s_max-(.9/n_components)),
                          (s_max,s_min+(.9/n_components)+(line_offset)*(n_components-1)/n_components), color="k", linewidth=1, linestyle=":")
        fig.add_artist(line)
    
    cax = fig.add_subplot(gs[j,j])
    cax.set_visible(False)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap), ax=cax, aspect = 10, fraction = 0.8)
    cbar.ax.tick_params(labelsize=fontsize)
    if use_metadata_coloring:
        cbar.ax.set_ylabel(column_color, fontsize=fontsize_label, labelpad=10)
    
    else:
        if coloring_property == "Missingness":
            cbar.ax.set_ylabel("# not NaN", fontsize=fontsize_label, labelpad=10)
        elif coloring_property == "Mean_Intensities":
            cbar.ax.set_ylabel("logarithmic\nmean intensities", fontsize=fontsize_label, labelpad=10)
        elif coloring_property == "Median_Intensities":
            cbar.ax.set_ylabel("logarithmic\nmedian intensities", fontsize=fontsize_label, labelpad=10)
        elif coloring_property=="Sum_Intensities":
            cbar.ax.set_ylabel("logarithmic\nsum of intesities", fontsize=fontsize_label, labelpad=10)
    if not continuous:
        cbar.ax.set_yticks(list(values_dict.values()))
        cbar.ax.set_yticklabels(values_dict.keys())
    
    plt.savefig(outfile, bbox_inches="tight")
    if return_fitted_data:
        if data_processed is not None:
            return data_origin_transformed, data_processed_transformed, list(data_origin_pca.columns), list(data_processed_pca.columns)
        else:
            return data_origin_transformed, data_origin_pca.columns
    return

def plot_pca_metadata_vectors(metadata, data_origin_transformed, data_origin_columns, data_processed_transformed = None, data_processed_columns = None, n_components = 4, min_samples = 3, max_entries = 10, max_entries_forced = 30, categorical_properties = None, verbose = False, top_n = None, noise_level = None, size = (10,10), fontsize = 15, fontsize_title = 20, line_offset = 0.025, fontsize_label = 12, outfile = "PCA_metadata_vectors.png"):
    """
    Plots vector of each metadata property with respect to its contribution driving the variance of the first n_components principal components.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties) with the same indices as the columns of data_origin.
    data_origin_transformed : np.ndarray
        Sample coordinates of the principal components of data_origin (Indices: Samples, Columns: Components).
    data_origin_columns : list
        Column names of data_origin, in the row order of data_origin_transformed.
    data_processed_transformed : np.ndarray, optional
        Sample coordinates of the principal components of data_processed (Indices: Samples, Columns: Components). The default is None.
    data_processed_columns : list, optional
        Column names of data_processed, in the row order of data_processed_transformed. The default is None.
    n_components : int, optional
        First (highest Eigenvalue) n_components the properties are projected on. The default is 4.
    min_samples : int, optional
        Minimum number of samples with a finite value in both the property and the component, below which the coordinate is set to NaN. The default is 3.
    max_entries : int, optional
        Maximum number of levels (excluding missing values) of a categorical property, above which the property is dropped. The default is 10.
    max_entries_forced : int, optional
        Maximum number of levels of a property matched by categorical_properties, above which the property is dropped. Set generously to catch free text. The default is 30.
    categorical_properties : tuple of str, optional
        Substrings, matched case-insensitively against the column names, marking properties that are one-hot encoded regardless of their dtype and limited by max_entries_forced instead of max_entries. An exact column name works as well, since it is a substring of itself. The default is None.
    verbose : bool, optional
        If True, the dropped properties and the numeric properties encoded as categories are printed. The default is False.
    top_n : int, optional
        If set, only the top_n vectors with the highest length are drawn within each pair of components. The default is None.
    noise_level : float, optional
        If set, vectors with a length below noise_level are not drawn. The length is the fraction of the variance of a property explained by the two components, so a value of 0.3 hides everything below 9% explained variance. The default is None.
    size : tuple of float, optional
        Size of the plot. The default is (10,10).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    fontsize_title : int, optional
        Fontsize of the title. The default is 15.
    line_offset : float, optional
        Offset of diagonal line (needs to adjusted if a lot of components are plotted. The default is 0.025.
    fontsize_label : int, optional
        Font size of the vector labels. The default is 12.
    outfile : str, optional
        File path for the figure. The default is "PCA_metadata_vectors.png".
    """

    vectors_origin = get_metadata_vectors(metadata, data_origin_transformed, data_origin_columns, n_components = n_components, min_samples = min_samples, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, verbose = verbose)
    vectors_processed = None
    if (data_processed_transformed is not None) and (data_processed_columns is not None):
        vectors_processed = get_metadata_vectors(metadata, data_processed_transformed, data_processed_columns, n_components = n_components, min_samples = min_samples, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, verbose = verbose)

    n_components = min(n_components, vectors_origin.shape[1])
    #properties, labels, colors = get_metadata_property_colors(vectors_origin.index, metadata)
    visible = []

    # Plotting
    fig = plt.figure(figsize=(n_components*2.5,n_components*2.5))
    gs = GridSpec(n_components,n_components, hspace= 0.1, wspace = 0.1)

    for i in range(n_components-1):
        for j in range(i+1,n_components):
            ax = fig.add_subplot(gs[j,i])
            entries = draw_vectors(ax, vectors_origin, metadata, i, j, fontsize_label, top_n = top_n, noise_level = noise_level)
            visible += [entry for entry in entries if entry not in visible]

            ax.tick_params(labelsize = fontsize)
            if i != 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("PC "+str(j+1), fontsize = fontsize, labelpad = 5)
            if j != n_components-1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("PC "+str(i+1), fontsize = fontsize, labelpad = 5)

            ax.set_axisbelow(True)
            ax.grid(alpha=0.5)

            if vectors_processed is not None:
                ax = fig.add_subplot(gs[i,j])
                entries = draw_vectors(ax, vectors_processed, metadata, j, i, fontsize_label, top_n = top_n, noise_level = noise_level)
                visible += [entry for entry in entries if entry not in visible]

                ax.tick_params(labelsize = fontsize, labeltop=True, labelright=True, labelbottom=False, labelleft = False)
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position("top")

                if j != n_components-1:
                    ax.set_yticklabels([])
                else:
                    ax.set_ylabel("PC "+str(i+1), fontsize = fontsize, labelpad = 5)
                if i != 0:
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel("PC "+str(j+1), fontsize = fontsize, labelpad = 5)

                ax.set_axisbelow(True)
                ax.grid(alpha=0.5)

    if vectors_processed is not None:
        s = 2.5*(np.log2(n_components)+1)/n_components
        s_min = 0.01*n_components
        s_max = 1 - s_min
        fig.text(1.5/(s*n_components),(s*n_components-1.5)/(s*n_components), "Before", rotation=315, va="center", ha="right",color="royalblue", fontsize = fontsize_title)
        fig.text(1.5/(s*n_components),(s*n_components-1.5)/(s*n_components), "After", rotation=315, va="center", ha="left",color="orange", fontsize = fontsize_title)
        line = plt.Line2D((s_min,s_max-(.9/n_components)),
                          (s_max,s_min+(.9/n_components)+(line_offset)*(n_components-1)/n_components), color="k", linewidth=1, linestyle=":")
        fig.add_artist(line)

    #lax = fig.add_subplot(gs[n_components-1,n_components-1])
    #lax.set_axis_off()
    #drawn = [name for entry, name in zip(vectors_origin.index, properties) if entry in visible]
    #handles = [plt.Line2D([0],[0], color=colors[name], linewidth=2, label=name) for name in dict.fromkeys(drawn)]
    #lax.legend(handles=handles, loc="center", fontsize=fontsize_label, frameon=False)

    plt.savefig(outfile, bbox_inches="tight")
    return


def plot_pca_metadata_eta2(metadata, data_origin_transformed, data_origin_columns, data_processed_transformed = None, data_processed_columns = None, n_components = 4, min_samples = 3, max_entries = 10, max_entries_forced = 30, categorical_properties = None, adjusted = False, verbose = False, top_n = None, noise_level = None, merge_duplicates = False, set_axis_lim = False, axis_offset = 0.1, size = (10,10), fontsize = 15, fontsize_title = 20, line_offset = 0.025, fontsize_label = 12, outfile = "PCA_metadata_eta2.png",):
    """
    Plots the squared correlation ratio (eta2) of each metadata property with respect to its contribution driving the variance of the first n_components principal components.
 
    All levels of a property are summarised by a single value in the range [0,1], so a property is a point within the unit square (relationship square) instead of a vector within the unit circle. 
    The distance from the origin measures how much of the plotted plane a property accounts for, the position along an axis how the effect is distributed over the two components.
 
    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties) with the same indices as the columns of data_origin.
    data_origin_transformed : np.ndarray
        Sample coordinates of the principal components of data_origin (Indices: Samples, Columns: Components).
    data_origin_columns : list
        Column names of data_origin, in the row order of data_origin_transformed.
    data_processed_transformed : np.ndarray, optional
        Sample coordinates of the principal components of data_processed (Indices: Samples, Columns: Components). The default is None.
    data_processed_columns : list, optional
        Column names of data_processed, in the row order of data_processed_transformed. The default is None.
    n_components : int, optional
        First (highest Eigenvalue) n_components the properties are projected on. The default is 4.
    min_samples : int, optional
        Minimum number of samples with a finite value in both the property and the component, below which the value is set to NaN. The default is 3.
    max_entries : int, optional
        Maximum number of levels (excluding missing values) of a categorical property, above which the property is dropped. The default is 10.
    max_entries_forced : int, optional
        Maximum number of levels of a property matched by categorical_properties, above which the property is dropped. Set generously to catch free text. The default is 30.
    categorical_properties : tuple of str, optional
        Substrings, matched case-insensitively against the column names, marking properties that are one-hot encoded regardless of their dtype and limited by max_entries_forced instead of max_entries. An exact column name works as well, since it is a substring of itself. The default is None.
    adjusted : bool, optional
        If True, eta2 is corrected for the number of levels of a property, which makes properties with many levels comparable to properties with few levels and can become slightly negative. The default is False.
    verbose : bool, optional
        If True, the dropped properties and the numeric properties encoded as categories are printed. The default is False.
    top_n : int, optional
        If set, only the top_n properties with the highest distance from the origin are drawn within each pair of components. The default is None.
    noise_level : float, optional
        If set, properties with a distance from the origin below noise_level are not drawn. The default is None.
    merge_duplicates : bool, optional
        If True, properties with identical values are merged into a single entry whose name lists all of them, since a property that is contained twice in the metadata (a batch as a name and as an identifier) explains exactly the same variance and would be drawn on top of itself. The default is False.
    set_axis_lim : bool, optional
        If True, the axes are limited to the range of the values instead of the full unit square, which separates the properties if all of them explain only a small fraction of the variance. One limit is calculated per component, so that it stays constant along the whole row and column of a dataset, and all components share the same span, so that the panels stay square. The default is False.
    axis_offset : float, optional
        Fraction of the span added on both sides if set_axis_lim == True. The default is 0.1.
    size : tuple of float, optional
        Size of the plot. The default is (10,10).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    fontsize_title : int, optional
        Fontsize of the title. The default is 15.
    line_offset : float, optional
        Offset of diagonal line (needs to adjusted if a lot of components are plotted. The default is 0.025.
    fontsize_label : int, optional
        Font size of the point labels. The default is 12.
    outfile : str, optional
        File path for the figure. The default is "PCA_metadata_eta2.png".
    """
 
    eta2_origin = get_metadata_eta2(metadata, data_origin_transformed, data_origin_columns, n_components = n_components, min_samples = min_samples, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, adjusted = adjusted, verbose = verbose)
    eta2_processed = None
    if (data_processed_transformed is not None) and (data_processed_columns is not None):
        eta2_processed = get_metadata_eta2(metadata, data_processed_transformed, data_processed_columns, n_components = n_components, min_samples = min_samples, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, adjusted = adjusted, verbose = verbose)
 
    if merge_duplicates:
        eta2_origin = merge_duplicate_entries(eta2_origin)
        if eta2_processed is not None:
            eta2_processed = merge_duplicate_entries(eta2_processed)
 
    n_components = min(n_components, eta2_origin.shape[1])
    lower = min(0.0, float(np.nanmin(eta2_origin.to_numpy())) if eta2_origin.notna().any().any() else 0.0)
    limits_origin = get_axis_limits(eta2_origin, n_components, axis_offset = axis_offset) if set_axis_lim else None
    limits_processed = None
    if eta2_processed is not None:
        limits_processed = get_axis_limits(eta2_processed, n_components, axis_offset = axis_offset) if set_axis_lim else None
 
    # Plotting
    fig = plt.figure(figsize=(n_components*2.5,n_components*2.5))
    gs = GridSpec(n_components,n_components, hspace= 0.1, wspace = 0.1)
 
    for i in range(n_components-1):
        for j in range(i+1,n_components):
            ax = fig.add_subplot(gs[j,i])
            draw_eta2(ax, eta2_origin, metadata, i, j, fontsize_label, top_n = top_n, noise_level = noise_level, lower = lower, limits = limits_origin)
 
            ax.tick_params(labelsize = fontsize)
            if i != 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("PC "+str(j+1), fontsize = fontsize, labelpad = 5)
            if j != n_components-1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("PC "+str(i+1), fontsize = fontsize, labelpad = 5)
 
            ax.set_axisbelow(True)
            ax.grid(alpha=0.5)
 
            if eta2_processed is not None:
                ax = fig.add_subplot(gs[i,j])
                draw_eta2(ax, eta2_processed, metadata, j, i, fontsize_label, top_n = top_n, noise_level = noise_level, lower = lower, limits = limits_processed)
 
                ax.tick_params(labelsize = fontsize, labeltop=True, labelright=True, labelbottom=False, labelleft = False)
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position("top")
 
                if j != n_components-1:
                    ax.set_yticklabels([])
                else:
                    ax.set_ylabel("PC "+str(i+1), fontsize = fontsize, labelpad = 5)
                if i != 0:
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel("PC "+str(j+1), fontsize = fontsize, labelpad = 5)
 
                ax.set_axisbelow(True)
                ax.grid(alpha=0.5)
 
    if eta2_processed is not None:
        s = 2.5*(np.log2(n_components)+1)/n_components
        s_min = 0.01*n_components
        s_max = 1 - s_min
        fig.text(1.5/(s*n_components),(s*n_components-1.5)/(s*n_components), "Before", rotation=315, va="center", ha="right",color="royalblue", fontsize = fontsize_title)
        fig.text(1.5/(s*n_components),(s*n_components-1.5)/(s*n_components), "After", rotation=315, va="center", ha="left",color="orange", fontsize = fontsize_title)
        line = plt.Line2D((s_min,s_max-(.9/n_components)),
                          (s_max,s_min+(.9/n_components)+(line_offset)*(n_components-1)/n_components), color="k", linewidth=1, linestyle=":")
        fig.add_artist(line)
 
    plt.savefig(outfile, bbox_inches="tight")
    return

def get_metadata_vectors(metadata, data_origin_transformed, data_origin_columns, data_processed_transformed = None, data_processed_columns = None, n_components = 4, min_samples = 3, max_entries = 10, max_entries_forced = 30, categorical_properties = None, verbose = False):
    """
    Calculation of the vector of each metadata property with respect to its contribution driving the variance of the first n_components principal components.
    Each coordinate is the Pearson correlation between the property and the sample coordinates of a principal component for supplementary variables. 
    Non-numeric properties are one-hot encoded beforehand and every level is projected separately as "Property=Level". 
    Properties with two levels are represented by a single vector.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties) with the same indices as the columns of data_origin.
    data_origin_transformed : np.ndarray
        Sample coordinates of the principal components of data_origin (Indices: Samples, Columns: Components).
    data_origin_columns : list
        Column names of data_origin, in the row order of data_origin_transformed.
    data_processed_transformed : np.ndarray, optional
        Sample coordinates of the principal components of data_processed (Indices: Samples, Columns: Components). The default is None.
    data_processed_columns : list, optional
        Column names of data_processed, in the row order of data_processed_transformed. The default is None.
    n_components : int, optional
        First (highest Eigenvalue) n_components the properties are projected on. The default is 4.
    min_samples : int, optional
        Minimum number of samples with a finite value in both the property and the component, below which the coordinate is set to NaN. The default is 3.
    max_entries : int, optional
        Maximum number of levels (excluding missing values) of a categorical property, above which the property is dropped. The default is 10.
    max_entries_forced : int, optional
        Maximum number of levels of a property matched by categorical_properties, above which the property is dropped. Set generously to catch free text. The default is 30.
    categorical_properties : tuple of str, optional
        Substrings, matched case-insensitively against the column names, marking properties that are one-hot encoded regardless of their dtype and limited by max_entries_forced instead of max_entries. An exact column name works as well, since it is a substring of itself. The default is None.
    verbose : bool, optional
        If True, the dropped properties and the numeric properties encoded as categories are printed. The default is False.

    Returns
    -------
    vectors_origin : pd.DataFrame
        Coordinates of each property of data_origin (Indices: Properties/Levels, Columns: PC_1 ... PC_n) in the range [-1,1].
    vectors_processed : pd.DataFrame, optional
        Coordinates of each property of data_processed (Indices: Properties/Levels, Columns: PC_1 ... PC_n) in the range [-1,1]. Only returned if data_processed_transformed and data_processed_columns are provided.

    """
    vectors_origin = get_vectors(metadata, data_origin_transformed, data_origin_columns, n_components = n_components, min_samples = min_samples, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, verbose = verbose)
    if (data_processed_transformed is None) or (data_processed_columns is None):
        return vectors_origin
    vectors_processed = get_vectors(metadata, data_processed_transformed, data_processed_columns, n_components = n_components, min_samples = min_samples, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, verbose = verbose)
    return vectors_origin, vectors_processed

def get_metadata_eta2(metadata, data_origin_transformed, data_origin_columns, data_processed_transformed = None, data_processed_columns = None, n_components = 4, min_samples = 3, max_entries = 10, max_entries_forced = 30, categorical_properties = None, adjusted = False, verbose = False):
    """
    Calculation of the squared correlation ratio (eta2) of each metadata property with respect to its contribution driving the variance of the first n_components principal components.

    Eta2 is the fraction of the variance of a principal component that is explained by a property, i.e. the coefficient of determination of the component regressed on all levels of that property. 
    Unlike the correlation vectors, all levels of a property are summarised by a single value in the range [0,1], which makes properties of different types and with a different number of levels comparable.
    This allows to track the effect size of a property before and after a processing step.
    For a numeric property eta2 is the squared Pearson correlation, for a property with two levels it is the squared coordinate of get_metadata_vectors.

    Eta2 increases with the number of levels of a property even for random labels, since one parameter is fitted per level. 
    adjusted = True applies the usual correction for the degrees of freedom, which makes properties with many levels Batch) comparable to properties with few levels (Sex) and can become slightly negative if a property explains no more than noise.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties) with the same indices as the columns of data_origin.
    data_origin_transformed : np.ndarray
        Sample coordinates of the principal components of data_origin (Indices: Samples, Columns: Components).
    data_origin_columns : list
        Column names of data_origin, in the row order of data_origin_transformed.
    data_processed_transformed : np.ndarray, optional
        Sample coordinates of the principal components of data_processed (Indices: Samples, Columns: Components). The default is None.
    data_processed_columns : list, optional
        Column names of data_processed, in the row order of data_processed_transformed. The default is None.
    n_components : int, optional
        First (highest Eigenvalue) n_components the properties are projected on. The default is 4.
    min_samples : int, optional
        Minimum number of samples with a finite value in both the property and the component, below which the value is set to NaN. The default is 3.
    max_entries : int, optional
        Maximum number of levels (excluding missing values) of a categorical property, above which the property is dropped. The default is 10.
    max_entries_forced : int, optional
        Maximum number of levels of a property matched by categorical_properties, above which the property is dropped. Set generously, since the number of batches of a study is not known beforehand, but low enough to still catch free text. The default is 30.
    categorical_properties : tuple of str, optional
        Substrings, matched case-insensitively against the column names, marking properties that are one-hot encoded regardless of their dtype and limited by max_entries_forced instead of max_entries. An exact column name works as well, since it is a substring of itself. The default is None.
    adjusted : bool, optional
        If True, eta2 is corrected for the number of levels of a property. The default is False.
    verbose : bool, optional
        If True, the dropped properties and the numeric properties encoded as categories are printed. The default is False.

    Returns
    -------
    eta2_origin : pd.DataFrame
        Squared correlation ratio of each property of data_origin (Indices: Properties, Columns: PC_1 ... PC_n) in the range [0,1].
    eta2_processed : pd.DataFrame, optional
        Squared correlation ratio of each property of data_processed (Indices: Properties, Columns: PC_1 ... PC_n) in the range [0,1]. Only returned if data_processed_transformed and data_processed_columns are provided.

    """
    eta2_origin = get_eta2(metadata, data_origin_transformed, data_origin_columns, n_components = n_components, min_samples = min_samples, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, adjusted = adjusted, verbose = verbose)
    if (data_processed_transformed is None) or (data_processed_columns is None):
        return eta2_origin
    eta2_processed = get_eta2(metadata, data_processed_transformed, data_processed_columns, n_components = n_components, min_samples = min_samples, max_entries = max_entries, max_entries_forced = max_entries_forced, categorical_properties = categorical_properties, adjusted = adjusted, verbose = verbose)
    return eta2_origin, eta2_processed

def draw_vectors(ax, vectors, metadata, i, j, fontsize_label, top_n = None, noise_level = None):
    """
    Drawing of the correlation circle of two principal components into a single axis.

    Parameters
    ----------
    ax : obj
        Axis the vectors are drawn into.
    vectors : pd.DataFrame
        Coordinates of each property (Indices: Properties/Levels, Columns: PC_1 ... PC_n).
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties).
    i : int
        Component plotted along the x-axis.
    j : int
        Component plotted along the y-axis.
    fontsize_label : int
        Font size of the vector labels.
    top_n : int, optional
        If set, only the top_n vectors with the highest length are drawn. The default is None.
    noise_level : float, optional
        If set, vectors with a length below noise_level are not drawn. The default is None.

    Returns
    -------
    entries : list
        Rows of vectors that have been drawn.

    """
    properties, labels, colors = get_metadata_property_colors(vectors.index, metadata)
    entries = get_visible_entries(vectors, i, j, top_n = top_n, noise_level = noise_level)

    ax.set_xlim(-1.25,1.25)
    ax.set_ylim(-1.25,1.25)
    ax.set_xticks([-1,0,1])
    ax.set_yticks([-1,0,1])
    ax.set_xticks([-0.5,0.5], minor=True)
    ax.set_yticks([-0.5,0.5], minor=True)
    ax.tick_params(which="minor", length=0)
    ax.grid(True, alpha=0.5, which="minor", color=plt.rcParams["grid.color"], linestyle=plt.rcParams["grid.linestyle"], linewidth=plt.rcParams["grid.linewidth"])
    ax.set_aspect("equal")

    ax.add_patch(plt.Circle((0,0), 1, edgecolor="k", facecolor="None", linewidth=1))
    ax.axhline(0, color="k", linewidth=0.5)
    ax.axvline(0, color="k", linewidth=0.5)

    points = []
    texts = []
    label_colors = []
    for entry in entries:
        index = list(vectors.index).index(entry)
        property_name = properties[index]
        label = labels[index]
        x = vectors.loc[entry,"PC_%i" % (i+1)]
        y = vectors.loc[entry,"PC_%i" % (j+1)]
        ax.annotate("", xy=(x,y), xytext=(0,0), arrowprops=dict(arrowstyle="-|>", color=colors[property_name], linewidth=1.5, shrinkA=0, shrinkB=0))
        points.append((x,y))
        texts.append(property_name if label == property_name else "%s : %s" % (property_name, label))
        label_colors.append(colors[property_name])

    draw_repelled_labels(ax, points, texts, label_colors, fontsize_label)
    return entries

def draw_eta2(ax, eta2, metadata, i, j, fontsize_label, top_n = None, noise_level = None, lower = 0.0, limits = None):
    """
    Drawing of the relationship square of two principal components into a single axis.
 
    Parameters
    ----------
    ax : obj
        Axis the points are drawn into.
    eta2 : pd.DataFrame
        Squared correlation ratio of each property (Indices: Properties, Columns: PC_1 ... PC_n).
    metadata : pd.DataFrame
        Metadata for each sample (Indices: Samples, Columns: Properties).
    i : int
        Component plotted along the x-axis.
    j : int
        Component plotted along the y-axis.
    fontsize_label : int
        Font size of the point labels.
    top_n : int, optional
        If set, only the top_n properties with the highest distance from the origin are drawn. The default is None.
    noise_level : float, optional
        If set, properties with a distance from the origin below noise_level are not drawn. The default is None.
    lower : float, optional
        Lower limit of both axes, below 0 for adjusted values. Ignored if limits is given. The default is 0.0.
    limits : list of tuple of float, optional
        Lower and upper limit of each component as returned by get_axis_limits. If None, the full unit square is shown. The default is None.
 
    Returns
    -------
    entries : list
        Rows of eta2 that have been drawn.
 
    """
    properties, labels, colors = get_metadata_property_colors(eta2.index, metadata)
    entries = get_visible_entries(eta2, i, j, top_n = top_n, noise_level = noise_level)
 
    if limits is None:
        ax.set_xlim(lower-0.05,1.05)
        ax.set_ylim(lower-0.05,1.05)
        ax.set_xticks([0,1])
        ax.set_yticks([0,1])
        ax.set_xticks([0.5], minor=True)
        ax.set_yticks([0.5], minor=True)
    else:
        ax.set_xlim(limits[i])
        ax.set_ylim(limits[j])
        ax.xaxis.set_major_locator(MaxNLocator(nbins = 3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins = 3))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="minor", length=0)
    ax.grid(True, alpha=0.5, which="minor", color=plt.rcParams["grid.color"], linestyle=plt.rcParams["grid.linestyle"], linewidth=plt.rcParams["grid.linewidth"])
    ax.set_aspect("equal")
 
    points = []
    texts = []
    label_colors = []
    for entry in entries:
        index = list(eta2.index).index(entry)
        property_name = properties[index]
        label = labels[index]
        x = eta2.loc[entry,"PC_%i" % (i+1)]
        y = eta2.loc[entry,"PC_%i" % (j+1)]
        ax.scatter(x, y, color=colors[property_name], s=60, zorder=3)
        points.append((x,y))
        texts.append(property_name if label == property_name else "%s : %s" % (property_name, label))
        label_colors.append(colors[property_name])
 
    draw_repelled_labels(ax, points, texts, label_colors, fontsize_label)
    return entries

def draw_repelled_labels(ax, points, texts, label_colors, fontsize_label, padding = 1.5, marker_size = 6.0):
    """
    Drawing of labels next to their anchor, avoiding the other labels and the markers.

    The size of each label is estimated from its number of characters and the font size and converted into data coordinates.
    Every label is then placed on the first free position out of eight directions around its anchor and three distances, which is the classical placement of a point label and needs no iteration.
    The labels are placed in the order they are given, so the most important one should come first, and a thin line is drawn as soon as a label could only be placed further away.

    Parameters
    ----------
    ax : obj
        Axis the labels are drawn into. The limits and the aspect have to be set beforehand.
    points : list of tuple of float
        Anchor of each label in data coordinates.
    texts : list of str
        Text of each label.
    label_colors : list
        Color of each label.
    fontsize_label : int
        Font size of the labels.
    padding : float, optional
        Additional space around each label in points. The default is 1.5.
    marker_size : float, optional
        Radius of a marker in points, kept free around each anchor. The default is 6.0.

    Returns
    -------
    None

    """
    if len(points) == 0:
        return

    origin = ax.transData.inverted().transform((0.0,0.0))
    unit = ax.transData.inverted().transform((1.0,1.0))
    pixel = ax.get_figure().dpi/72.0
    per_x = abs(unit[0]-origin[0])*pixel
    per_y = abs(unit[1]-origin[1])*pixel

    widths = np.array([(0.55*len(text)*fontsize_label + 2*padding)*per_x for text in texts])
    heights = np.array([(1.15*fontsize_label + 2*padding)*per_y for text in texts])

    anchors = np.array(points, dtype = float)
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    marker_x = marker_size*per_x
    marker_y = marker_size*per_y

    directions = [(0,1),(1,0),(0,-1),(-1,0),(1,1),(1,-1),(-1,1),(-1,-1)]
    distances = [1.0,1.8,2.6]

    positions = np.zeros_like(anchors)
    scales = np.ones(len(anchors))
    for index in range(len(anchors)):
        best_position = None
        best_cost = None
        best_distance = 1.0
        for distance in distances:
            for direction_x, direction_y in directions:
                position = np.array([anchors[index,0] + direction_x*(widths[index]/2.0+marker_x)*distance,
                                     anchors[index,1] + direction_y*(heights[index]/2.0+marker_y)*distance])
                cost = 0.0
                cost += max(0.0, (x_min+widths[index]/2.0)-position[0]) + max(0.0, position[0]-(x_max-widths[index]/2.0))
                cost += max(0.0, (y_min+heights[index]/2.0)-position[1]) + max(0.0, position[1]-(y_max-heights[index]/2.0))
                cost = cost*1000.0
                for other in range(index):
                    overlap_x = (widths[index]+widths[other])/2.0 - abs(position[0]-positions[other,0])
                    overlap_y = (heights[index]+heights[other])/2.0 - abs(position[1]-positions[other,1])
                    if (overlap_x > 0) and (overlap_y > 0):
                        cost += overlap_x*overlap_y/(widths[index]*heights[index])
                for other in range(len(anchors)):
                    overlap_x = widths[index]/2.0 + marker_x - abs(position[0]-anchors[other,0])
                    overlap_y = heights[index]/2.0 + marker_y - abs(position[1]-anchors[other,1])
                    if (overlap_x > 0) and (overlap_y > 0):
                        cost += overlap_x*overlap_y/(widths[index]*heights[index])
                if (best_cost is None) or (cost < best_cost):
                    best_cost = cost
                    best_position = position
                    best_distance = distance
            if best_cost < 1e-9:
                break
        positions[index] = best_position
        scales[index] = best_distance

    for index, text in enumerate(texts):
        if scales[index] > 1.0:
            ax.plot([anchors[index,0],positions[index,0]], [anchors[index,1],positions[index,1]], color=label_colors[index], linewidth=0.5, alpha=0.6, zorder=2)
        ax.text(positions[index,0], positions[index,1], text, color=label_colors[index], fontsize=fontsize_label, ha="center", va="center", zorder=4)
    return

def plot_distribution(data_origin, metadata, data_processed=None, sample_type = None, sample_type_column = None, step="", outfile = "Distribution.png", size = (10,5), nbins = 50, smoothing_factor = 1e-2, fontsize = 15, xlabel = "log2(PQ)", cmap = "Spectral"):
    """
    Plotting of the intensity distribution per sample.

    Parameters
    ----------
    data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame
	    Metadata for each sample with the same indices as the columns of data_origin.
	data_processed : pd.DataFrame, optional
	    Data (Indices: Precursors, Columns: Samples) after the preprocessing step.
    sample_type : str, optional
	    Sample type in Metadata[sample_type_column] for which indices should be extracted. The default is None.
    sample_type_column : str, optional
	    Column in Metadata that distinguishes sample type. If None, samples will not be distinguished. The default is None.
    step : str, optional
        Preprocessing step (used for the title). The default is "".
    outfile : str, optional
        File path for the figure. The default is "Distribution.png".
    size : tuple of float, optional
        Size of the plot. The default is (10,5).
    nbins : int, optional
        Number of bins for estimating the distribution. The default is 50.
    smoothing_factor : str, optional
        Smoothing factor for interpolating the distribution based on a histogram. The default is 1e-2.
    fontsize : fontsize : int, optional
        Font size of the labels. The default is 15.
    xlabel : str, optional
        x-axis label. The default is "log2(PQ)".
    cmap : str, optional
        Colormap used for coloring. The default is "Spectral".
    """
    
    fig = plt.figure(figsize=size)
    fig.set_size_inches(size)
    
    if data_processed is not None:
        gs = GridSpec(1,2, width_ratios = [0.5,0.5], hspace = 0.05, wspace = 0.05)
    else:
        gs = GridSpec(1,1, hspace= 0.05)
    
    ax_origin = fig.add_subplot(gs[0,0])
    
    if data_processed is not None:
        max_val = max(data_origin.max().max(),data_processed.max().max())
        min_val = min(data_origin.min().min(),data_processed.min().min())
    else:
        max_val = data_origin.max().max()
        min_val = data_origin.min().min()
        
    bin_edges = np.linspace(min_val, max_val, nbins + 1)
    width = bin_edges[1] - bin_edges[0]
    
    if sample_type_column:
        if not sample_type:
            raise ValueError("sample_type is not defined.")
        indices_origin = [index for index in metadata.index if metadata.loc[index,sample_type_column]==sample_type and index in data_origin.columns]
        if data_processed is not None:
            indices_processed = [index for index in metadata.index if metadata.loc[index,sample_type_column]==sample_type and index in data_processed.columns]
    else:
        indices_origin = [index for index in metadata.index if index in data_origin.columns]
        if data_processed is not None:
            indices_processed = [index for index in metadata.index if index in data_processed.columns]
    
    #distribution_origin = [np.histogram(data_origin[column],bins=nbins,range=(min_val,max_val),density=True)[0] for column in indices_origin]
    
    X = data_origin.loc[:, indices_origin].to_numpy()
    bin_idx = np.full(X.shape, -1, dtype=np.int64)
    valid = np.isfinite(X)
    bin_idx[valid] = ((X[valid] - min_val) / width).astype(np.int64)
    mask = (bin_idx >= 0) & (bin_idx < nbins)

    ncols = X.shape[1]
    col_idx = np.broadcast_to(np.arange(ncols),X.shape)
    flat_bins = bin_idx[mask] + col_idx[mask] * nbins
    counts = np.bincount(flat_bins, minlength=ncols * nbins).reshape(ncols, nbins)
    densities = counts / (X.shape[0] * width)
    distribution_origin = densities.tolist()
    
    xs = np.linspace(min_val+(max_val-min_val)/(2*nbins),max_val-(max_val-min_val)/(2*nbins),nbins)
    xs_fine = np.linspace(min_val+(max_val-min_val)/(2*nbins),max_val-(max_val-min_val)/(2*nbins),250)
    
    cmap_ = plt.get_cmap(cmap)
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    color_values_origin = np.linspace(0,1,len(indices_origin))
    colors_origin = [cmap_(norm(el)) for el in color_values_origin] 
    
    max_y = 0
    
    for item, color in zip(distribution_origin,colors_origin):
        spl = UnivariateSpline(xs,item)
        spl.set_smoothing_factor(smoothing_factor)
        y_mod = spl(xs_fine)
        ax_origin.plot(xs_fine,y_mod,color=color,alpha=0.5)
        if max(y_mod) > max_y:
            max_y = max(y_mod)
            
    ax_origin.tick_params(labelsize = fontsize)
    ax_origin.set_xlabel(xlabel,fontsize=fontsize)
    ax_origin.set_ylabel("Density",fontsize=fontsize)
    ax_origin.grid(axis="y")
    ax_origin.set_xlim(min_val,max_val)
    ax_origin.set_ylim(0,max_y*1.05)
    
    if data_processed is not None:
        ax_processed = fig.add_subplot(gs[0,1])
        
        #distribution_processed = [np.histogram(data_processed[column],bins=nbins,range=(min_val,max_val),density=True)[0] for column in indices_processed]
        X = data_processed.loc[:, indices_processed].to_numpy()
        bin_idx = np.full(X.shape, -1, dtype=np.int64)
        valid = np.isfinite(X)
        bin_idx[valid] = ((X[valid] - min_val) / width).astype(np.int64)
        mask = (bin_idx >= 0) & (bin_idx < nbins)

        ncols = X.shape[1]
        col_idx = np.broadcast_to(np.arange(ncols),X.shape)
        flat_bins = bin_idx[mask] + col_idx[mask] * nbins
        counts = np.bincount(flat_bins, minlength=ncols * nbins).reshape(ncols, nbins)
        densities = counts / (X.shape[0] * width)
        distribution_processed = densities.tolist()
        
        color_values_processed = np.linspace(0,1,len(indices_processed))
        colors_processed = [cmap_(norm(el)) for el in color_values_processed] 
        for item, color in zip(distribution_processed,colors_processed):
            spl = UnivariateSpline(xs,item)
            spl.set_smoothing_factor(smoothing_factor)
            y_mod = spl(xs_fine)
            ax_processed.plot(xs_fine,y_mod,color=color,alpha=0.5)
            if max(y_mod) > max_y:
                max_y = max(y_mod)   
        ax_processed.set_xlim(min_val,max_val)
        ax_processed.set_ylim(0,max_y*1.05)
        
        ax_origin.set_title("Before "+step, fontsize=fontsize)
        ax_processed.set_title("After "+step, fontsize=fontsize)
    
        ax_processed.tick_params(labelsize = fontsize)
        ax_processed.set_xlabel(xlabel,fontsize=fontsize)
        ax_processed.set_yticklabels([])
        ax_processed.grid(axis="y")
    
    plt.savefig(outfile, bbox_inches="tight")
    return

### TODO: Add that also a list of peptides can be provided instead of plot_by
def plot_scatter_series(data_origin, metadata, column_metadata, data_processed = None, sample_type = None, sample_type_column = None, plot_by = "Correlation", outfile = "Scatter_series.png", data_type = "continuous", top_n = 3, size = (5,3), ylabel = "log2(Prec. Abundance)", fontsize = 15, use_metadata_coloring = False, cmap = None, continuous = False, column_color = None, plot_linreg = False, base = 2, step="", min_samples=10, limit_y_label=None):
    """
    Plotting of Intensity versus covariate for top_n precursors of plot_by property.

    Parameters
    ----------
    data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame
	    Metadata for each sample with the same indices as the columns of data_origin.
    column_metadata : str
        Column in Metadata containing the covariate.
    data_processed : pd.DataFrame, optional
	    Data (Indices: Precursors, Columns: Samples) after the preprocessing step.
    sample_type : str, optional
	    Sample type in Metadata[sample_type_column] for which indices should be extracted. The default is None.
    sample_type_column : str, optional
	    Column in Metadata that distinguishes sample type. If None, samples will not be distinguished. The default is None.
    plot_by : str, optional
	       Property for which data should be calculated. Options: "Missingness", "Mean", "Median", "CV", "norm_MAE", "norm_MSE", "Correlation", "Corr_pvalue". The default is "Correlation".
    outfile : str, optional
        File path for the figure. The default is "Scatter_series.png".
    data_type : str, optional
	    If plot_by == "Correlation" or plot_by == "Corr_pvalue", type of data that is correlated. Options: "continuous", "ordinal", "nominal". The default is "continuous".
    top_n : int, optional
        top_n IDs of plot_by to be plotted. The default is 3.
    size : tuple of float, optional
        Size of the plot. The default is (5,3).
    ylabel : str, optional
        y-axis label. The default is "log2(Prec. Abundance)".
    fontsize : int, optional
        Font size of the labels. The default is 15.
    use_metadata_coloring : bool, optional
        If True, Metadata are used for coloring. The default is False.
    cmap : str, optional
        Colormap used for coloring. The default is None.
    continuous : bool, optional
        If True continuous data are assumed for coloring e.g. Order. If False discrete data are assumed e.g. Batch. The default is False.
    column_color : str, optional
        If use_metadata_coloring == True, Metadata column used for coloring. The default is None.
    plot_linreg : bool, optional
        If True linear regression over all intensities of selected peptides is plotted. The default is False.
    base : float, optional
	    Base of the logarithmic transformation. The default is 2.
    step : str, optional
        Preprocessing step (used for the title). The default is "".
    min_samples : int, optional
        Minimal numbers of samples to plot. The default is 10.
    limit_y_label : int,optional
        Maximal number of characters for the ylabel.
    """
        
    if column_metadata:
        try:
            metadata[column_metadata] = pd.to_numeric(metadata[column_metadata])
        except:
            raise TypeError("Meta data column needs to be numeric.")
    else: 
        metadata = metadata.copy()
        column_metadata = "Order"
        metadata[column_metadata] = np.arange(0,metadata.shape[0])
    
    try:
        metadata[column_metadata] = pd.to_numeric(metadata[column_metadata])
    except:
        raise TypeError("Meta data column needs to be numeric.")
    
    # sample_type selection
    if data_processed is not None:
        indices_origin, indices_processed = get_indices_by_type(data_origin, metadata, data_processed, sample_type_column, sample_type)
    else:
        indices_origin = get_indices_by_type(data_origin, metadata, data_processed, sample_type_column, sample_type)
    
    # Calculate property
    property_peptide = get_property(data_origin, indices_origin, metadata, column_metadata=column_metadata, plot_by=plot_by, data_type=data_type, take_abs=True, min_samples=min_samples)
    # Get n top hits
    selection = property_peptide.sort_values(ascending=False).index[:top_n]
    
    if data_processed is not None:
        size = (size[0]*2,(size[1]*top_n)+1)
    else:
        size = (size[0],(size[1]*top_n)+1)
    
    if not cmap:
        if continuous:
            cmap = "viridis"
        else:
            cmap = "Spectral"
    
    if use_metadata_coloring:
        if data_processed is not None:
            colors_origin, colors_processed, norm, values_dict = get_colors_by_metadata(metadata, column_color, indices_origin, indices_processed=indices_processed, cmap=cmap, continuous=continuous)
        else:
            colors_origin, norm, values_dict = get_colors_by_metadata(metadata, column_color, indices_origin, indices_processed=None, cmap=cmap, continuous=continuous)

    fig = plt.figure(figsize=size)
    if data_processed is not None:
        gs = GridSpec(top_n,2, width_ratios = [0.5,0.5], wspace = 0.05, hspace = 0.05)
    else:
        gs = GridSpec(top_n,1, hspace = 0.05)
    
    x_lim = (np.min(metadata[column_metadata].values),np.max(metadata[column_metadata].values))
    
    for ind, el in enumerate(selection):
        
        # estimate y_lim
        if data_processed is not None:
            y_lim = (np.min([np.min(data_origin.loc[el,indices_origin]),np.min(data_processed.loc[el,indices_processed])])*0.98,
                     np.max([np.max(data_origin.loc[el,indices_origin]),np.max(data_processed.loc[el,indices_processed])])*1.02)
        else:
            y_lim = (np.min(data_origin.loc[el,indices_origin])*0.98,np.max(data_origin.loc[el,indices_origin])*1.02)
            
        ax_scatter_origin = fig.add_subplot(gs[ind,0])
        if use_metadata_coloring:
            ax_scatter_origin.scatter(metadata.loc[indices_origin,column_metadata], data_origin.loc[el,indices_origin], c=colors_origin)
        else:
            ax_scatter_origin.scatter(metadata.loc[indices_origin,column_metadata], data_origin.loc[el,indices_origin])
        
        if not ind == top_n-1:
            ax_scatter_origin.set_xticklabels([])
        else:
            ax_scatter_origin.set_xlabel(column_metadata, fontsize=fontsize)
        ax_scatter_origin.set_xlim(x_lim)
        ax_scatter_origin.set_ylim(y_lim)
        ax_scatter_origin.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        if limit_y_label:
            if len(el) > limit_y_label:
                el_ = el[:limit_y_label]+"..."
            else:
                el_ = el
        else:
            el_ = el
        ax_scatter_origin.set_ylabel(ylabel+"\n"+el_, fontsize=fontsize)
        ax_scatter_origin.set_axisbelow(True)
        ax_scatter_origin.grid(axis="y")
        ax_scatter_origin.tick_params(labelsize=fontsize)
        
        if plot_linreg:
            if use_metadata_coloring:
                for value in set(metadata.loc[indices_origin,column_color]):
                    indices_origin_nonan = [index for index in indices_origin if not np.isnan(metadata.loc[index,column_metadata]) and not np.isnan(data_origin.loc[el,index]) and metadata.loc[index,column_color]==value]
                    x = [np.min(metadata.loc[indices_origin_nonan,column_metadata]),np.max(metadata.loc[indices_origin_nonan,column_metadata])]
                    coef = np.polyfit(metadata.loc[indices_origin_nonan,column_metadata], data_origin.loc[el,indices_origin_nonan],1) #### NaN Handling
                    linreg = np.poly1d(coef)
                    ax_scatter_origin.plot(x,linreg(x),c="k")
            else:
                indices_origin_nonan = [index for index in indices_origin if not np.isnan(metadata.loc[index,column_metadata]) and not np.isnan(data_origin.loc[el,index])]
                x = x_lim
                coef = np.polyfit(metadata.loc[indices_origin_nonan,column_metadata], data_origin.loc[el,indices_origin_nonan],1) #### NaN Handling
                linreg = np.poly1d(coef)
                ax_scatter_origin.plot(x,linreg(x),c="k")
        
        if data_processed is not None:
            ax_scatter_processed = fig.add_subplot(gs[ind,1])
            if use_metadata_coloring:
                ax_scatter_processed.scatter(metadata.loc[indices_processed,column_metadata], data_processed.loc[el,indices_processed], c=colors_processed)
            else:
                ax_scatter_processed.scatter(metadata.loc[indices_processed,column_metadata], data_processed.loc[el,indices_processed])
            
            if not ind == top_n-1:
                ax_scatter_processed.set_xticklabels([])
            else:
                ax_scatter_processed.set_xlabel(column_metadata, fontsize=fontsize)
            ax_scatter_processed.set_xlim(x_lim)    
            ax_scatter_processed.set_ylim(y_lim)
            ax_scatter_processed.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax_scatter_processed.set_yticklabels([])
            ax_scatter_processed.yaxis.tick_right()
            ax_scatter_processed.set_axisbelow(True)
            ax_scatter_processed.grid(axis="y")
            ax_scatter_processed.tick_params(labelsize=fontsize)
            
            if ind == 0:
                ax_scatter_origin.set_title("Before "+step, fontsize=fontsize)
                ax_scatter_processed.set_title("After "+step, fontsize=fontsize)
            
            if plot_linreg:
                if use_metadata_coloring:
                    for value in set(metadata.loc[indices_processed,column_color]):
                        indices_processed_nonan = [index for index in indices_processed if not np.isnan(metadata.loc[index,column_metadata]) and not np.isnan(data_processed.loc[el,index]) and metadata.loc[index,column_color]==value]
                        x = [np.min(metadata.loc[indices_processed_nonan,column_metadata]),np.max(metadata.loc[indices_processed_nonan,column_metadata])]
                        coef = np.polyfit(metadata.loc[indices_processed_nonan,column_metadata], data_processed.loc[el,indices_processed_nonan],1) #### NaN Handling
                        linreg = np.poly1d(coef)
                        ax_scatter_processed.plot(x,linreg(x),c="k")
                else:
                    indices_processed_nonan = [index for index in indices_processed if not np.isnan(metadata.loc[index,column_metadata]) and not np.isnan(data_processed.loc[el,index])]
                    x = x_lim
                    coef = np.polyfit(metadata.loc[indices_processed_nonan,column_metadata], data_processed.loc[el,indices_processed_nonan],1) #### NaN Handling
                    linreg = np.poly1d(coef)
                    ax_scatter_processed.plot(x,linreg(x),c="k")
    
    plt.savefig(outfile, bbox_inches="tight")
    return

def plot_decision_boundary_estimation_curve(data_origin, metadata, method="response-binned-splines", outfile="Decision_boundary_estimation.png", sample_type_column=None, sample_type=None, size=(5,5), n_knots=3, boundary=0.5, detection_proportion_cap=1, x_lim=None, x_label="Mean(log2(Abundance))", plot_all_fits=True, fontsize=15):
    """
    Plotting of the decision boundary detection curve for missingness type estimation.

    Parameters
    ----------
    data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame
	    Metadata for each sample with the same indices as the columns of data_origin.
    method : str, optional
        method used for fitting. Options: "logit-linear", "logit-linear-fit", "logit-linear-splines", "logit-linear-capped", "logit-linear-biascorrected", "response-binned-splines". The default is "response-binned-splines".
    outfile : str, optional
        File path for the figure. The default is "Decision_boundary_estimation.png".
    sample_type_column : str, optional
	    Column in Metadata that distinguishes sample type. If None, samples will not be distinguished. The default is None.
    sample_type : str, optional
	    Sample type in Metadata[sample_type_column] for which indices should be extracted. The default is None.
    size : tuple of float, optional
        Size of the plot. The default is (5,5).
    n_knots : int, optional
        Number of knots for spline calculation. The default is 3.
    boundary : float, optional
        Completeness boundary for decision making. The default is 0.5.
    detection_proportion_cap : float, optional
        Scaling for logit-linear-capped. The default is 1.
    x_lim : tuple of float, optional
        x-axis limit. The default is None.
    x_label : str, optional
        x-axis label. The default is "Mean(log2(Abundance))".
    plot_all_fits : bool, optional
        If True all fits are plotted. If False only the fit for the selected method is plotted. The default is True.
    fontsize : int, optional
        Font size of the labels. The default is 15.
    """
    
    def func1_logitlinear(x,b0,b1):
            return b0 + b1 * x
        
    def func2_logitlinear_df(X,b):
        return np.sum(X*b,axis=1)

    def func3_capped_logitlinear(x,b0,b1):
        return detection_proportion_cap*np.exp(b0+b1*x)/(1+np.exp(b0+b1*x))

    def func4_logitlinear_biascorrected(X,b0,b1):
        return b0 + b1*X['Avg Log2 Abundance'] - (0.5*(b1**2)*X['Var Log2 Abundance'])
    
    indices_origin = get_indices_by_type(data_origin, metadata, sample_type_column=sample_type_column, sample_type = sample_type)
        
    #compute means, variances and missingness
    means = data_origin[indices_origin].mean(axis=1) 
    variances = data_origin[indices_origin].var(axis=1,ddof=1) 
    missingness = data_origin[indices_origin].isnull().mean(axis=1)
    
    #dataframe containing results
    res_df = pd.DataFrame({
        'Avg Log2 Abundance': means,
        'Var Log2 Abundance': variances,
        'Completeness': 1 - missingness,
        'Logit Completeness': np.log((1 - missingness)/missingness)
    }).replace([-np.inf,np.inf],np.nan).dropna(subset=["Avg Log2 Abundance","Completeness"])
    #sort res_df by abundance
    res_df.sort_values(by="Avg Log2 Abundance",inplace=True)
    
    cmap = mpl.colormaps['viridis']
    
    N=6
    
    fits = {"Completeness pred_logit-linear":cmap(0/(N-1)),
            "Completeness pred_logit-linear-fit":cmap(1/(N-1)),
            "Completeness pred_logit-linear-splines":cmap(2/(N-1)),
            "Completeness pred_logit-linear-capped":cmap(3/(N-1)),
            "Completeness pred_logit-linear-biascorrected":cmap(4/(N-1)), 
            "Completeness pred_response-binned-splines":cmap(5/(N-1))}
    
    #### Curve fits ####
    #Linear model with binomial distribution for logit completeness
    X = sm.add_constant(res_df['Avg Log2 Abundance'])
    model = sm.GLM(res_df['Completeness'], X, family=sm.families.Binomial())
    result = model.fit()
    intercept, coef = result.params
    res_df["Completeness pred_logit-linear"] = 1/(1+np.exp(-func1_logitlinear(res_df["Avg Log2 Abundance"],intercept,coef)))
    
    popt, pcov = curve_fit(func1_logitlinear,res_df.dropna(subset='Logit Completeness')['Avg Log2 Abundance'],res_df.dropna(subset='Logit Completeness')['Logit Completeness'],method="trf",loss="cauchy")
    res_df["Completeness pred_logit-linear-fit"] = 1/(1+np.exp(-func1_logitlinear(res_df['Avg Log2 Abundance'],*popt)))
    
    indices = res_df.dropna(subset=["Avg Log2 Abundance","Completeness","Logit Completeness"]).index
    X = sm.add_constant(res_df.loc[indices,'Avg Log2 Abundance'])
    spline = SplineTransformer(n_knots=n_knots, degree=1)
    splines = spline.fit_transform(res_df.loc[indices,'Avg Log2 Abundance'].to_numpy().reshape(-1,1),res_df.loc[indices,"Logit Completeness"])
    X = sm.add_constant(splines)
    model = sm.GLM(res_df.loc[indices,'Completeness'], X, family=sm.families.Binomial())
    result = model.fit()
    popt = result.params
    res_df.loc[indices,"Completeness pred_logit-linear-splines"] = 1/(1+np.exp(-func2_logitlinear_df(X=X,b=popt.to_numpy().reshape(1,-1))))
    
    popt, pcov = curve_fit(func3_capped_logitlinear,res_df['Avg Log2 Abundance'],res_df['Completeness'],bounds=([-np.inf,0],[np.inf,1]))
    res_df["Completeness pred_logit-linear-capped"] = func3_capped_logitlinear(res_df['Avg Log2 Abundance'],*popt)
    
    indices = res_df.dropna().index
    popt, pcov = curve_fit(func4_logitlinear_biascorrected,res_df.loc[indices],res_df.loc[indices,'Logit Completeness'])
    logit_completeness = func4_logitlinear_biascorrected(res_df.dropna(),*popt)
    popt, pcov = curve_fit(func1_logitlinear,res_df.loc[indices,"Avg Log2 Abundance"],logit_completeness)
    res_df.loc[indices,"Completeness pred_logit-linear-biascorrected"] = 1/(1+np.exp(-func1_logitlinear(res_df.loc[indices,"Avg Log2 Abundance"],*popt)))
    
    # Spline the response variable (Detection Probability) using pd.qcut()
    res_df['response_bin'] = pd.qcut(res_df['Completeness'], q=100, duplicates='drop')
    # Compute mean abundance and mean detection probability per bin
    grouped = res_df.groupby('response_bin',observed=False).agg(abundance_mean=('Avg Log2 Abundance', 'median'),detection_mean=('Completeness', 'median')).dropna()
    # Prepare data for GLM
    X = sm.add_constant(grouped['abundance_mean'])
    y = grouped['detection_mean']
    # Fit binomial GLM (using proportions directly with freq_weights)
    model = sm.GLM(y, X, family=sm.families.Binomial())
    result = model.fit()
    # Decision boundary calculation
    intercept, coef = result.params
    res_df["Completeness pred_response-binned-splines"] = 1/(1+np.exp(-func1_logitlinear(res_df["Avg Log2 Abundance"],intercept,coef)))
    #### ####
    
    decision_boundary = res_df.loc[(res_df["Completeness pred_"+method].dropna()>=boundary).idxmax(),"Avg Log2 Abundance"]
    
    fig,ax = plt.subplots()
    fig.set_size_inches(size)
    
    ax.scatter(means,1 - missingness,alpha=1,c="w",zorder=1)
    ax.scatter(means,1 - missingness,alpha=0.1,zorder=2)
    
    if not x_lim:
        x_lim = (np.min(means)-0.2,np.max(means)+0.2)
    
    for fit, color in fits.items():
        if fit.split("_")[-1] == method:
            ax.plot(res_df["Avg Log2 Abundance"],res_df[fit],ls="-",alpha=1,lw=2, c=color, label=fit.split("_")[-1],zorder=20)
        elif plot_all_fits:
            ax.plot(res_df["Avg Log2 Abundance"],res_df[fit],ls="--",alpha=1,lw=2, c=color, label=fit.split("_")[-1],zorder=3)
    
    ax.plot(x_lim,[boundary,boundary],c="k",ls=":",zorder=4)
    ax.plot([decision_boundary,decision_boundary],[0,1],c="k",ls=":",zorder=5)
        
    ax.set_xlim(x_lim)
    ax.set_ylim(0,1)
    ax.tick_params(labelsize=fontsize)
    
    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylabel("Completeness", fontsize=fontsize)
    ax.set_title("Decision Boundary: {:.2f}".format(decision_boundary), fontsize=fontsize)
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fontsize)
    ax.set_axisbelow(True)
    ax.grid()
    
    plt.savefig(outfile,bbox_inches="tight")
    
def plot_correlation_distribution(data_origin, metadata, data_processed = None, sample_type = None, sample_type_column = None, column_metadata = "Order", data_type = "continuous", seed = 42, nbins = 50, size = (10,5), fontsize = 15, fontsize_legend = 15, outfile = "Correlation_distribution.png", use_same_random = True, legend_loc = "upper right", step = ""):
    """
    Plotting of correlation distribution of covariate versus intensities compared to a random distribution (randomized covariate).

    Parameters
    ----------
    data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame
	    Metadata for each sample with the same indices as the columns of data_origin.
	data_processed : pd.DataFrame, optional
	    Data (Indices: Precursors, Columns: Samples) after the preprocessing step.
    sample_type : str, optional
	    Sample type in Metadata[sample_type_column] for which indices should be extracted. The default is None.
    sample_type_column : str, optional
	    Column in Metadata that distinguishes sample type. If None, samples will not be distinguished. The default is None.
    column_metadata : str, optional
	    Column in Metadata containing the covariate for correlation. The default is "Order".
    data_type : data_type : str, optional
	    If plot_by == "Correlation" or plot_by == "Corr_pvalue", type of data that is correlated. Options: "continuous", "ordinal", "nominal". The default is "continuous".
    seed : int, optional
	    Seed for shuffeling the data. The default is 42.
    nbins : int, optional
        Number of bins for the histogram. The default is 50.
    size : tuple of float, optional
        Size of the plot. The default is (10,5).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    fontsize_legend : int, optional
        Font size of the legend. The default is 15.
    outfile : str, optional
        File path for the figure. The default is "Correlation_distribution.png".
    use_same_random : bool, optional
        If True the same randomized covariate is used for comparison before and after preprocessing step. The default is True.
    legend_loc : str, optional
        Location of the legend. The default is "upper right".
    step : str, optional
        Preprocessing step (used for the title).. The default is "".
    """
    
    ### change to dummy conversion in future????
    try:
        metadata[column_metadata] = pd.to_numeric(metadata[column_metadata])
    except:
        raise TypeError("Meta data column needs to be numeric.")
    
    # sample_type selection
    if data_processed is not None:
        indices_origin, indices_processed = get_indices_by_type(data_origin, metadata, data_processed, sample_type_column, sample_type)
    else:
        indices_origin = get_indices_by_type(data_origin, metadata, data_processed, sample_type_column, sample_type)
        
    # Calculate correlation
    corr_peptide = get_property(data_origin, indices_origin, metadata, column_metadata=column_metadata, data_type=data_type, take_abs=False)
    corr_peptide_random = get_property(data_origin, indices_origin, metadata, column_metadata=column_metadata, data_type=data_type, shuffle=True, seed=seed)
    
    range_min_origin = np.min((np.min(corr_peptide),np.min(corr_peptide_random)))
    range_max_origin = np.max((np.max(corr_peptide),np.max(corr_peptide_random)))
    
    hist_corr_origin = np.histogram(corr_peptide,bins=nbins,range=(range_min_origin,range_max_origin))
    hist_corr_random_origin = np.histogram(corr_peptide_random,bins=nbins,range=(range_min_origin,range_max_origin))
    
    if data_processed is not None:
        # Calculate correlation
        corr_peptide = get_property(data_processed, indices_processed, metadata, column_metadata=column_metadata, data_type=data_type, take_abs=False)
        if not use_same_random:
            corr_peptide_random = get_property(data_processed, indices_processed, metadata, column_metadata=column_metadata, data_type=data_type, shuffle=True, seed=seed)
        
        range_min_processed = np.min((np.min(corr_peptide),np.min(corr_peptide_random)))
        range_max_processed = np.max((np.max(corr_peptide),np.max(corr_peptide_random)))
    
        hist_corr_processed = np.histogram(corr_peptide,bins=nbins,range=(range_min_origin,range_max_origin))
        if not use_same_random:
            hist_corr_random_processed = np.histogram(corr_peptide_random,bins=nbins,range=(range_min_origin,range_max_origin))
        else: 
            hist_corr_random_processed = hist_corr_random_origin
            
    if data_processed is not None:
        x_min = np.min((range_min_origin,range_min_processed))
        x_max = np.max((range_max_origin,range_max_processed))
        y_max = np.max((np.max((np.max(hist_corr_origin[0]),np.max(hist_corr_random_origin[0]))),np.max((np.max(hist_corr_processed[0]),np.max(hist_corr_random_processed[0])))))
    else:
        x_min = range_min_origin
        x_max = range_max_origin
        y_max = np.max((np.max(hist_corr_origin[0]),np.max(hist_corr_random_origin[0])))
        
    fig = plt.figure(figsize=size)
    fig.set_size_inches(size)
    
    if data_processed is not None:
        gs = GridSpec(1,2, width_ratios = [0.5,0.5], hspace = 0.05, wspace = 0.05)
    else:
        gs = GridSpec(1,1, hspace= 0.05)
        
    ax_origin = fig.add_subplot(gs[0,0])
    ax_origin.bar((hist_corr_origin[1][:-1]+hist_corr_origin[1][1:])/2, hist_corr_origin[0],alpha=1, color="w",width=hist_corr_origin[1][1]-hist_corr_origin[1][0])
    ax_origin.bar((hist_corr_random_origin[1][:-1]+hist_corr_random_origin[1][1:])/2, hist_corr_random_origin[0],alpha=1, color="w",width=hist_corr_random_origin[1][1]-hist_corr_random_origin[1][0])
    ax_origin.bar((hist_corr_origin[1][:-1]+hist_corr_origin[1][1:])/2, hist_corr_origin[0],alpha=0.5, color="C0",width=hist_corr_origin[1][1]-hist_corr_origin[1][0],label=column_metadata)
    ax_origin.bar((hist_corr_random_origin[1][:-1]+hist_corr_random_origin[1][1:])/2, hist_corr_random_origin[0],alpha=0.5, color="C1",width=hist_corr_random_origin[1][1]-hist_corr_random_origin[1][0],label="Random")
    
    ax_origin.tick_params(labelsize = fontsize)
    if data_type == "continuous":
        ax_origin.set_xlabel("Pearson R",fontsize=fontsize)
    elif data_type == "ordinal":
        ax_origin.set_xlabel("Spearman R",fontsize=fontsize)
    elif data_type == "nominal":
        ax_origin.set_xlabel("Point biserial R",fontsize=fontsize)
    ax_origin.set_ylabel("Occurence",fontsize=fontsize)
    ax_origin.set_axisbelow(True)
    ax_origin.grid(axis="y")
    ax_origin.set_xlim(x_min,x_max)
    ax_origin.set_ylim(0,y_max*1.05)
    
    if data_processed is not None:
        ax_processed = fig.add_subplot(gs[0,1])
        ax_processed.bar((hist_corr_processed[1][:-1]+hist_corr_processed[1][1:])/2, hist_corr_processed[0],alpha=1, color="w",width=hist_corr_processed[1][1]-hist_corr_processed[1][0])
        ax_processed.bar((hist_corr_random_processed[1][:-1]+hist_corr_random_processed[1][1:])/2, hist_corr_random_processed[0],alpha=1, color="w",width=hist_corr_random_processed[1][1]-hist_corr_random_processed[1][0])
        ax_processed.bar((hist_corr_processed[1][:-1]+hist_corr_processed[1][1:])/2, hist_corr_processed[0],alpha=0.5, color="C0",width=hist_corr_processed[1][1]-hist_corr_processed[1][0])
        ax_processed.bar((hist_corr_random_processed[1][:-1]+hist_corr_random_processed[1][1:])/2, hist_corr_random_processed[0],alpha=0.5, color="C1",width=hist_corr_random_processed[1][1]-hist_corr_random_processed[1][0])
        
        ax_processed.tick_params(labelsize = fontsize)
        if data_type == "continuous":
            ax_processed.set_xlabel("Pearson R",fontsize=fontsize)
        elif data_type == "ordinal":
            ax_processed.set_xlabel("Spearman R",fontsize=fontsize)
        elif data_type == "nominal":
            ax_processed.set_xlabel("Point biserial R",fontsize=fontsize)
        ax_processed.set_ylabel("")
        ax_processed.set_yticklabels([])
        ax_processed.set_axisbelow(True)
        ax_processed.grid(axis="y")
        ax_processed.set_xlim(x_min,x_max)
        ax_processed.set_ylim(0,y_max*1.05)
        
        ax_origin.set_title("Before "+step, fontsize=fontsize)
        ax_processed.set_title("After "+step, fontsize=fontsize)
        
    try:
        ax_origin.legend(loc=legend_loc, fontsize=fontsize_legend)
    except:
        ax_origin.legend(loc="best", fontsize=fontsize_legend)
        print("legend_loc not found. Using best")
    
    plt.savefig(outfile, bbox_inches="tight")
    return

def plot_imputation_vs_covariate(data_origin, metadata, mask, column_metadata = "Order", sample_type = None, sample_type_column = None, size = (5,5), fontsize = 15, fontsize_legend = 15, outfile = "Imputation.png", detailed = True, imputation_value = 99):
    """
    Plotting of imputation versus covariate.

    Parameters
    ----------
    data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame
	    Metadata for each sample with the same indices as the columns of data_origin.
    mask : pd.DataFrame
        Imputation mask (detailed = True)/Filter mask (detailed = False).
    column_metadata : str, optional
	    Column in Metadata containing the covariate. The default is "Order".
    sample_type : str, optional
	    Sample type in Metadata[sample_type_column] for which indices should be extracted. The default is None.
    sample_type_column : str, optional
	    Column in Metadata that distinguishes sample type. If None, samples will not be distinguished. The default is None.
    size : tuple of float, optional
        Size of the plot. The default is (5,5).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    fontsize_legend : int, optional
        Font size of the legend. The default is 15.
    outfile : str, optional
        File path for the figure. The default is "Imputation.png".
    detailed : bool, optional
        If True plot each missingness type separately. The default is True.
    imputation_value : int, optional
        Identifier in filter mask for imputet values. The default is 99.
    """
    
    try:
        metadata[column_metadata] = pd.to_numeric(metadata[column_metadata])
    except:
        raise TypeError("Meta data column needs to be numeric.")
    
    # sample_type selection
    indices_origin = get_indices_by_type(data_origin, metadata, sample_type_column=sample_type_column, sample_type=sample_type)
    
    fig,ax = plt.subplots()
    fig.set_size_inches(size)
    x = metadata.loc[indices_origin,column_metadata]
    
    if detailed:
        mask.fillna("0",inplace=True)
        mask_entries = np.unique(mask.values)
        for entry in mask_entries:
            if entry!= "0":
                count_entry = (mask==entry).sum()
                ax.plot(x,count_entry[indices_origin],label=entry)
                
        imputation_comp = (mask!="0").sum()
        ax.plot(x,imputation_comp[indices_origin], c="k",ls="--",label="Total")
    else:
        if any(mask=="HA") or any(mask=="LA") or any(mask=="ProtAbsence") or any("ProtAbsence/HA") or any("ProtAbsence/LA"):
            mask = mask.isna().astype(float)*0.99
        temp_mask=mask.map(lambda x: int(str(float(x)).split(".")[-1]))  
        imputed_precs=(temp_mask!=imputation_value).sum(axis=0)
        ax.plot(x,imputed_precs[indices_origin], c="k",ls="--",label="Total")
    
    ax.tick_params(which="both",labelsize = fontsize)
    ax.set_xlabel(column_metadata, fontsize=fontsize)
    ax.set_ylabel("Occurence", fontsize=fontsize)
    if detailed:
        ax.set_ylim(0,np.max(imputation_comp[indices_origin])*1.1)
    else:
        ax.set_ylim(0,np.max(imputed_precs[indices_origin])*1.1)
    
    ax.legend(loc="center left",fontsize=fontsize_legend, bbox_to_anchor=(1, 0.5))
    
    ax.set_axisbelow(True)
    ax.grid()
    
    plt.savefig(outfile, bbox_inches="tight")
    return

def plot_distribution_hist(data, outfile = "Normalisation_factors.png", size = (5,5), fontsize = 15, xlabel = "Normalisation factors", transform = False, base = 2, nbins = 20, use_log_axis = False, baseline = None):
    """
    Plotting of distribution histogram for e.g. Normalization factors.

    Parameters
    ----------
    data : pd.Series
        Data to be plotted e.g. Normalization factors.
    outfile : str, optional
        File path for the figure. The default is "Normalisation_factors.png".
    size : tuple of float, optional
        Size of the plot. The default is (5,5).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    xlabel : str, optional
        x-axis label. The default is "Normalisation factors".
    transform : bool, optional
        If True, data are delogarithmize. The default is False.
    base : float, optional
	    Base of the logarithmic transformation. The default is 2.
    nbins : int, optional
        Number of bins for the histogram calculation. The default is 20.
    use_log_axis : bool, optional
        If True, y-axis scale is set to logarithmic. The default is False.
    baseline : int, optional
        Baseline of the distribution. The default is None.
    """
    
    if transform:
        data = base**data
    
    fig,ax = plt.subplots()
    fig.set_size_inches(size)
    
    data_hist = np.histogram(data,bins=nbins)
    ax.bar((data_hist[1][:-1]+data_hist[1][1:])/2, data_hist[0],alpha=1, color="C0",width=data_hist[1][1]-data_hist[1][0])
    
    if baseline is not None:
        ax.plot([baseline,baseline],[0,np.max(data_hist[0])+1],c="k",ls="--",lw=2)
    
    if use_log_axis:
        ax.set_yscale('log')
    
    ax.tick_params(axis="both", labelsize=fontsize)
    ax.set_ylabel("Occurence", fontsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    
    ax.set_axisbelow(True)
    ax.grid()
    
    plt.savefig(outfile,bbox_inches="tight")
    return

def plot_completeness_vs_abundance_scatter(data_origin, metadata, sample_type = None, sample_type_column = None, outfile = "Completness_vs_Abundance.png", size = (7.5,7.5), fontsize = 15, axis = 1, gridsize = 25, n_bins = 25, cmap = "inferno", bins_setup='log', x_label = "Mean(log2(Abundance))", y_label = "Completeness"):
    """
    Plotting of completeness versus abundance grouped as a hexagonal binned plot.

    Parameters
    ----------
    data_origin : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame
	    Metadata for each sample with the same indices as the columns of data_origin.
    sample_type : str, optional
	    Sample type in Metadata[sample_type_column] for which indices should be extracted. The default is None.
    sample_type_column : str, optional
	    Column in Metadata that distinguishes sample type. If None, samples will not be distinguished. The default is None.
    outfile : str, optional
        File path for the figure. The default is "Completness_vs_Abundance.png".
    size : tuple of float, optional
        Size of the plot. The default is (7.5,7.5).
    fontsize : int, optional
        Font size of the labels. The default is 15.
    axis : int, optional
        Axis for mean value calculations. The default is 1.
    gridsize : int, optional
        Number of the grid points per dimension. The default is 25.
    n_bins : int, optional
        Number of bins for histogram calculation. The default is 25.
    cmap : str, optional
        Colormap used for coloring. The default is "inferno".
    bins_setup : str or int, optional
        Discretization of binned values. Check matplotlib.axes.Axes.hexbin for detailed description. The default is 'log'.
    x_label : str, optional
        x-axis label. The default is "Mean(log2(Abundance))".
    y_label : str, optional
        y-axis label. The default is "Completeness".
    """
    
    # sample_type selection
    indices_origin = get_indices_by_type(data_origin, metadata, sample_type_column=sample_type_column, sample_type=sample_type)
    
    cmap_ = plt.get_cmap(cmap)
    color_hist = cmap_(0.5)
    
    fig = plt.figure()
    fig.set_size_inches(size)
    gs = GridSpec(3,4,height_ratios=[11/6,1/6,10], width_ratios=[10,2/3,2/3,2/3], hspace=0.01, wspace=0.01)
    
    abundance = data_origin[indices_origin].mean(axis=axis)
    completeness = 1-data_origin[indices_origin].isna().sum(axis=axis)/data_origin[indices_origin].shape[1]
    
    x_lim = (np.min(abundance),np.max(abundance))
    y_lim = (np.min(completeness),np.max(completeness))
    
    ax_hex = fig.add_subplot(gs[2,0])
    im = ax_hex.hexbin(abundance, completeness, gridsize=gridsize, cmap=cmap_, bins=bins_setup)
    
    ax_hex.tick_params(axis="both",labelsize=fontsize)
    ax_hex.set_xlabel(x_label, fontsize=fontsize)
    ax_hex.set_ylabel(y_label, fontsize=fontsize)
    ax_hex.set_xlim(x_lim)
    ax_hex.set_ylim(y_lim)
    
    ax_hex.set_axisbelow(True)
    ax_hex.grid()
    
    ax_hist_x = fig.add_subplot(gs[:2,0])
    ax_hist_x.axis("off")
    
    hist_abundance = np.histogram(abundance, range=x_lim, bins=n_bins, density=True)
    
    ax_hist_x.bar((hist_abundance[1][1:]+hist_abundance[1][:-1])/2, hist_abundance[0], color=color_hist, width=hist_abundance[1][1]-hist_abundance[1][0])
    ax_hist_x.set_xlim(x_lim)
    ax_hist_x.set_ylim(0,np.max(hist_abundance[0]))
    
    ax_hist_x.grid()
    
    ax_hist_y = fig.add_subplot(gs[2,1:])
    ax_hist_y.axis("off")
    
    hist_completeness = np.histogram(completeness, range=y_lim, bins=n_bins, density=True)
    
    ax_hist_y.barh((hist_completeness[1][1:]+hist_completeness[1][:-1])/2, hist_completeness[0], color=color_hist, height=hist_completeness[1][1]-hist_completeness[1][0])
    ax_hist_y.set_xlim(0,np.max(hist_completeness[0]))
    ax_hist_y.set_ylim(y_lim)
    
    ax_cbar = fig.add_subplot(gs[0,2])
    
    plt.colorbar(im, orientation='vertical', cax=ax_cbar)
    ax_cbar.tick_params(axis="y",labelsize=fontsize)
    for tick in ax_cbar.yaxis.get_majorticklabels():
        tick.set_verticalalignment("center")
        
    plt.savefig(outfile,bbox_inches="tight")

def plot_overall_completeness(data, metadata=None, column_metadata=None, size = (5,5), outfile = "Overall_Completeness.png", fontsize = 15, fontsize_legend = 11, x_label = "Samples", y_label = "Sorted precurssors", color_present = "#6F15BD", color_abscent = "#BD9915", legend_loc = "upper left", plot_limits = True, limits = [1,0.75,0.5,0.33]):
    """
    Plotting of completeness per sample. IDs are sorted according to completeness.

    Parameters
    ----------
    data : pd.DataFrame
	    Data (Indices: Precursors, Columns: Samples) before the preprocessing step.
	metadata : pd.DataFrame, optional
	    Metadata for each sample with the same indices as the columns of data_origin. The default is None.
    column_metadata : str, optional
	    Column in Metadata containing the order. If set to None, sample order according to data is used. The default is None.
    size : tuple of float, optional
        Size of the plot. The default is (5,5).
    outfile : str, optional
        File path for the figure. The default is "Overall_Completeness.png".
    fontsize : int, optional
        Font size of the labels. The default is 15.
    fontsize_legend : int, optional
        Font size of the legend. The default is 11.
    x_label : str, optional
        x-axis label. The default is "Samples".
    y_label : str, optional
        y-axis label. The default is "Sorted precurssors".
    color_present : str, optional
        Color either default names or Hex color code for present IDs. The default is "#6F15BD".
    color_abscent : TYPE, optional
        Color either default names or Hex color code for abscent IDs. The default is "#BD9915".
    legend_loc : str, optional
        Location of the legend. The default is "upper left".
    plot_limits : bool, optional
        If True plot horizontal lines for the set limits. The default is True.
    limits : list of floats, optional
        Completeness limits to be highlighted. The default is [1,0.75,0.5,0.33].
    """
    
    if column_metadata:
        try:
            order = pd.to_numeric(metadata[column_metadata])
        except:
            raise TypeError("Meta data column needs to be numeric.")
    else: 
        order = np.arange(0,data.shape[1])
    

    data_ = data.copy()
    data_.fillna(0,inplace=True)
    data_[data_!=0] = 1
    
    sums_ = data_.T.sum().sort_values()[::-1]
    data_ = data_.loc[sums_.index]
    
    if plot_limits:
        limits_indices = []
        for limit in limits:
            try:
                limits_indices.append(np.where(data_.index==sums_[sums_>=data_.shape[1]*limit].index[-1])[0][0])
            except:
                limits_indices.append(0)
        #limits_indices = [np.where(data_.index==sums_[sums_>=data_.shape[1]*limit].index[-1])[0][0] for limit in limits]
    
    limit_label = "Limits : "
    for ind,limit in enumerate(limits):
        if ind < len(limits)-1:
            limit_label += "{:.2f}".format(limit)+", "
        else:
            limit_label += "{:.2f}".format(limit)
    
    cmap = mpl.colors.ListedColormap([color_abscent,color_present])
    
    fig, ax = plt.subplots()
    fig.set_size_inches(size)
    ax.matshow(data_.iloc[:,order],cmap=cmap)
    
    if plot_limits:
        legend_elements = [Patch(facecolor=color_present, edgecolor=color_present, label="Present"),
                           Patch(facecolor=color_abscent, edgecolor=color_abscent, label="Absent"),
                           Line2D([0], [0], color='k', linestyle="--", label=limit_label)]
    else:
        legend_elements = [Patch(facecolor=color_present, edgecolor=color_present, label="Present"),
                           Patch(facecolor=color_abscent, edgecolor=color_abscent, label="Absent")]
    
    ax.legend(loc = legend_loc, handles = legend_elements, fontsize=fontsize_legend)
    
    if plot_limits:
        for limit in limits_indices:
            ax.plot([0,data_.shape[1]],[limit,limit],c="k",ls="--")
    
    ax.set_xlim(0,data_.shape[1])
    #if column_metadata:
    #  ax.set_xlabel(column_metadata, fontsize=fontsize)
    #else:
    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize)
    ax.tick_params(axis="both", labelsize=fontsize)        
    ax.axes.get_xaxis().set_ticklabels([])
    ax.set_title("Overall Completeness: {:.2f}".format(sum(sums_)/(data_.shape[0]*data_.shape[1])), fontsize=fontsize)
    set_aspect_equal(ax)
    
    plt.savefig(outfile, bbox_inches="tight")
    return

#%%

### Maybe move to another code to separate plotting functions and import plotting functions?

def init_module():
    """
    Module initialization.

    Returns
    -------
    workdir : str
        path of the working directory.
    step : str
        Name of the step evaluated.
    meta_df_fp : str
        File path of the Metadata.
    config_eval_fp : str
        File path of the evaluation configuration file.
    """
    
    parser = argparse.ArgumentParser(description="Automated evaluation")
    parser.add_argument("-w","--workingdir", help="working directory where the preprocessing structure was build", required=True)
    parser.add_argument("-m","--metapath", help="Filepath to metadata file", required=True) 
    parser.add_argument("-s","--step", help="step that should be evaluated", required=False, default="All") 
    parser.add_argument("-e","--config_eval", help="Filepath to evaluation configuration file", required=True) 
    parser.add_argument("--docker", action="store_true")
    args = parser.parse_args()
    
    workdir = args.workingdir.replace("\\","/")
    step = args.step
    meta_df_fp = args.metapath.replace("\\","/")
    config_eval_fp = args.config_eval.replace("\\","/")
    
    return workdir, step, meta_df_fp, config_eval_fp, args.docker

def get_config_dict(config_fp, keys, step):
    """
    Tranform evaluation configuration file to disctionary

    Parameters
    ----------
    config_fp : str
        File path of the evaluation configuration file..
    keys : list of str
        Figure keys.
    step : str
        Step to be evaluated. If all steps should be evaluated at once set to "All".

    Returns
    -------
    dict_config_setup : dict
        Dictionary containing the configuration setup for every step and figure.
    steps : list of str
        Steps to be evaluated later.
    """
    
    dict_config_setup = {}
    steps = []
    with open(config_fp,"r") as file:
        for line in file.read().splitlines():
            if line:
                if line.startswith("%%%%%%%%%%"):
                    key = re.search(("%%%%%%%%%%(.+?)%%%%%%%%%%"),line).group(1)
                    dict_config_setup.update({key:{}})
                    counter = -1
                    if step == "All" and key != "00_GENERAL":
                        steps.append(key)
                    elif step == key:
                        steps.append(key)
                elif line == "%-%-%-%":
                    counter += 1
                    if key == "00_GENERAL":
                        key_2 = "Setup"
                    else:
                        key_2 = keys[counter]
                    dict_config_setup[key].update({key_2:{}})
                elif line != "%-%-%-%":
                    dict_config_setup[key][key_2].update({line.split("=")[0].strip():line.split("=")[1].strip()})
    return dict_config_setup, steps
    
#%%
if __name__ == "__main__":
    if len(sys.argv) <= 1:
        
        step = input("Please define the step (All is also an option): ")
        if step[3:] not in ["FILTERING_LONG","FILTERING_WIDE","NORMALISATION","IMPUTATION","CORRECTION","SUMMARISATION"] and not step == "All":
            print("Step not found. No step is evaluated evaluated")
            sys.exit(0)

        workdir = askdirectory(title="Please select your working directory", initialdir=os.path.abspath(os.getcwd()))
        meta_df_fp = askopenfilename(title="Please select the metadata file", initialdir=workdir)       
        config_fp = askopenfilename(title="Please select the evaluation config file", initialdir=workdir)
        docker_setup = False
    
    else:
        workdir, step, meta_df_fp, config_fp, docker_setup = init_module()
    
    keys = ["Data","F01","F02","F03","F04","F05","F06","F07","F08","F09","F10","F11"]    
    dict_config_setup, steps = get_config_dict(config_fp, keys, step)
                
    parameters_dict_general = helper.load_parameters(dict_config_setup["00_GENERAL"]["Setup"]["param_fp"])
    try:
        metadata_df, _ = helper.load_metadata(metadata_fp=meta_df_fp,
                                              metadata_sheet=parameters_dict_general["metadata_sheet"], 
                                              filename_col=parameters_dict_general["metadata_filename_column"], 
                                              encoding=parameters_dict_general["metadata_encoding"])
        metadata_df.set_index(parameters_dict_general["metadata_sample_column"], inplace=True, drop=False)

    except:
        raise ValueError("Metadata "+meta_df_fp+" could not be loaded. Skip evaluation.")
        sys.exit(1)
        
    if not steps:
        raise ValueError("Step could not be found")
        
    if not docker_setup:
        try:
            captions = pd.read_csv(os.path.join(os.path.split(os.path.abspath(__file__))[0],"Figures_Captions.tsv"), sep = "\t", index_col = 0)
        except:
            captions = pd.read_csv("./Figures_Captions.tsv", sep = "\t", index_col = 0)
    else:
        captions = pd.read_csv("Figures_Captions.tsv", sep = "\t", index_col = 0)
    
#%%
    for step in steps:
    #step = "01_FILTERING_LONG"
    #step = "05_CORRECTION"
        data_processed = helper.load_data(dict_config_setup[step]["Data"]["data_processed"])
        parameters_dict = helper.load_parameters(dict_config_setup[step]["Data"]["param_fp"])
        #if step == "01_FILTERING_LONG":
        if dict_config_setup[step]["Data"]["data_origin"] == "None":
            data_origin = data_processed
        else:
            data_origin = helper.load_data(dict_config_setup[step]["Data"]["data_origin"])
        mask_processed = helper.load_mask(dict_config_setup[step]["Data"]["data_processed"])
        mask_description_processed = helper.load_mask_description(dict_config_setup[step]["Data"]["data_processed"])
        
        log = open(os.path.join(re.search("(.+?)F01",dict_config_setup[step]["F01"]["outfile"]).group(1),"Figure_log.txt"),"w")
        
        #%%   
        key_figure = "F01"
        if dict_config_setup[step][key_figure]["plot_int_box"] == "1":
            if "CORRECTION" in step and "batch_cols_meta" in parameters_dict:
                use_metadata_coloring = True
                column_color = parameters_dict["batch_cols_meta"][0]
            else:
                use_metadata_coloring = False
                column_color = None
            if "FILTERING_LONG" in step:
                data_processed_ = None
            else:
                data_processed_ = data_processed
            if "SUMMARISATION" in step:
                ylabel_int = "Unique IDs"
            else:
                ylabel_int = "# Precursors"
                
            if parameters_dict_general["sample_type"] is not None:
                addition = " "+parameters_dict_general["sample_type_column"]+":"+parameters_dict_general["sample_type"]+"."
            else:
                addition = " All."
            
            try:
                plot_qa_box_int_and_number(data_origin, 
                                           metadata_df, 
                                           data_processed_, 
                                           sample_type = parameters_dict_general["sample_type"], 
                                           sample_type_column = parameters_dict_general["sample_type_column"],
                                           size = ast.literal_eval(dict_config_setup[step][key_figure]["size"]), 
                                           step = step, 
                                           outfile = dict_config_setup[step][key_figure]["outfile"], 
                                           use_metadata_coloring = use_metadata_coloring, 
                                           column_color = column_color,
                                           ylabel_int = ylabel_int
                                           )
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition=addition)
            except:
                print("Error in boxplot plotting. Skipped.")
            
        #%%
        key_figure = "F02"
        if dict_config_setup[step][key_figure]["plot_filter_mask"] == "1":  
            try:
                plot_filter_mask(mask_processed,
                                 mask_description_processed, 
                                 outfile = dict_config_setup[step][key_figure]["outfile"]
                                 )
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition="")
            except:
                print("Error in filter mask plotting. Skipped.")
    
        #%%
        key_figure = "F03"
        if dict_config_setup[step][key_figure]["plot_cv"] == "1":
            
            if parameters_dict_general["sample_type_column"] is None:
                groups_ = []
                addition = " None."           
            else:
                addition = " "+parameters_dict_general["sample_type_column"]+"."
                groups_ = parameters_dict_general["sample_type_column"]
            
            try:
                append_o, append_p = plot_CV_violin(data_origin, 
                                                    metadata_df, 
                                                    data_processed=data_processed, 
                                                    group_columns = [groups_], 
                                                    outfile = dict_config_setup[step][key_figure]["outfile"],
                                                    y_lim=(0,2)
                                                    )
                
                if append_o:
                    addition += " At least 1 precursor in original data exceeds a CV of 2."
                if append_p:
                    addition += " At least 1 precursor in processed data exceeds a CV of 2."
                
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition=addition)
            except:
                print("Error in CV plotting. Skipped.")
            
        #%%
        key_figure = "F04"
        if dict_config_setup[step][key_figure]["plot_pca"] == "1":
            if "CORRECTION" in step:
                if parameters_dict["correct_method"] in ["LIMMA_R","COMBAT_R","HARMONY_R","harmony","combat"]:
                    try:
                        dict_config_setup[step][key_figure]["column_color"] = parameters_dict["batch_cols_meta"][0]
                        addition = " "+parameters_dict["batch_cols_meta"][0]+"."
                    except:
                        print("batch_cols_meta not defined.")
                    dict_config_setup[step][key_figure]["continuous"] = "False"
                elif parameters_dict["correct_method"] in ["adjust_confounders","adjust_confounders_batchwise"]:
                    try:
                        dict_config_setup[step][key_figure]["column_color"] = parameters_dict["adjust_cols"][0]
                        addition = " "+parameters_dict["adjust_cols"][0]+"."
                        col_type = parameters_dict["adjust_cols_dtypes"][0]
                        if col_type in "continuous":
                            dict_config_setup[step][key_figure]["continuous"] = "True"
                        else:
                            dict_config_setup[step][key_figure]["continuous"] = "False"
                    except:
                        print("adjust_cols not defined.")
                dict_config_setup[step][key_figure]["coloring_property"] = None
            elif "SUMMARISATION" in step:
                dict_config_setup[step][key_figure]["column_color"] = parameters_dict_general["sample_type_column"]
                addition = " "+str(parameters_dict["sample_type_column"])+"."
                dict_config_setup[step][key_figure]["coloring_property"] = None
            else:
                dict_config_setup[step][key_figure]["column_color"] = None
                addition = " "+dict_config_setup[step][key_figure]["coloring_property"]+"."
                
            try:
                plot_PCA(data_origin, 
                         metadata_df, 
                         data_processed = data_processed, 
                         outfile = dict_config_setup[step][key_figure]["outfile"], 
                         continuous = ast.literal_eval(dict_config_setup[step][key_figure]["continuous"]),
                         column_color = dict_config_setup[step][key_figure]["column_color"], 
                         use_metadata_coloring = ast.literal_eval(dict_config_setup[step][key_figure]["use_metadata_coloring"]),
                         coloring_property = dict_config_setup[step][key_figure]["coloring_property"]
                         )
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition=addition)
            except:
                print("Error in PCA plotting. Skipped.")
    
        #%%
        key_figure = "F05"
        if dict_config_setup[step][key_figure]["plot_distribution"] == "1":
            try:
                plot_distribution(data_origin, 
                                  metadata_df, 
                                  data_processed = data_processed, 
                                  outfile = dict_config_setup[step][key_figure]["outfile"],
                                  step = step
                                  )
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition="")
            except:
                print("Error in distribution plotting. Skipped.")
            
        #%%
        key_figure = "F06"
        if dict_config_setup[step][key_figure]["plot_scatter_series"] == "1":
            if "NORMALISATION" in step:
                dict_config_setup[step][key_figure]["column_metadata"] = None
                dict_config_setup[step][key_figure]["data_type"] = "continuous"
                dict_config_setup[step][key_figure]["continuous"] = True
                dict_config_setup[step][key_figure]["use_metadata_coloring"] = False
                dict_config_setup[step][key_figure]["column_color"] = None
                dict_config_setup[step][key_figure]["plot_linreg"] = True
                addition_1 = "Inferred Order (continuous)."
                addition_2 = ""
            elif "CORRECTION" in step:
                try:
                    dict_config_setup[step][key_figure]["column_metadata"] = parameters_dict["adjust_cols"][0]
                    dict_config_setup[step][key_figure]["data_type"] = parameters_dict["adjust_cols_dtypes"][0]
                    if parameters_dict["correct_method"] in ["adjust_confounders_batchwise"]:
                        dict_config_setup[step][key_figure]["use_metadata_coloring"] = True
                        dict_config_setup[step][key_figure]["column_color"] = parameters_dict["batch_cols_meta"][0]
                        dict_config_setup[step][key_figure]["plot_linreg"] = False
                        dict_config_setup[step][key_figure]["continuous"] = False
                        addition_2 = "Colored according to: "+parameters_dict["batch_cols_meta"][0]+"."
                    else:
                        dict_config_setup[step][key_figure]["use_metadata_coloring"] = False
                        dict_config_setup[step][key_figure]["column_color"] = None
                        dict_config_setup[step][key_figure]["plot_linreg"] = True
                        dict_config_setup[step][key_figure]["continuous"] = False
                        addition_2 = ""
                    addition_1 = " "+parameters_dict["adjust_cols"][0]+" ("+parameters_dict["adjust_cols_dtypes"][0]+")"+". "
                    #if dict_config_setup[step][key_figure]["data_type"] in "continuous":
                    #    dict_config_setup[step][key_figure]["continuous"] = True
                    #else:
                    #    dict_config_setup[step][key_figure]["continuous"] = False
                except:
                    print("adjust_cols not defined.")
            
            try:
                plot_scatter_series(data_origin, 
                                    metadata_df, 
                                    column_metadata = dict_config_setup[step][key_figure]["column_metadata"], 
                                    data_processed = data_processed, 
                                    data_type = dict_config_setup[step][key_figure]["data_type"], 
                                    use_metadata_coloring = dict_config_setup[step][key_figure]["use_metadata_coloring"], 
                                    continuous = dict_config_setup[step][key_figure]["continuous"], 
                                    column_color = dict_config_setup[step][key_figure]["column_color"],
                                    plot_linreg = dict_config_setup[step][key_figure]["plot_linreg"], 
                                    outfile = dict_config_setup[step][key_figure]["outfile"], 
                                    plot_by = dict_config_setup[step][key_figure]["plot_by"], 
                                    step = step,
                                    min_samples = 10,
                                    limit_y_label=15
                                    )
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition=addition_1+addition_2)
            except:
                print("Error in scatter series plotting. Skipped.")
        
        #%%
        key_figure = "F07"
        if dict_config_setup[step][key_figure]["plot_decision_boundary"] == "1":
            try:
                plot_decision_boundary_estimation_curve(data_origin, 
                                                        metadata_df, 
                                                        outfile= dict_config_setup[step][key_figure]["outfile"],
                                                        method = parameters_dict["intensity_dependent_missingness_fit_method"],
                                                        boundary = float(parameters_dict["intensity_dependent_missingness_decision_boundary"])
                                                        )
                addition = " "+parameters_dict["intensity_dependent_missingness_fit_method"]+"."
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition=addition)
            except:
                print("Error in decision boundary plotting. Skipped.")
            
        #%%
        key_figure = "F08"   
        if dict_config_setup[step][key_figure]["plot_correlation"] == "1":
            
            try:
                plot_correlation_distribution(data_origin, 
                                              metadata_df, 
                                              data_processed = data_processed, 
                                              column_metadata = parameters_dict["adjust_cols"][0], 
                                              data_type = parameters_dict["adjust_cols_dtypes"][0],  
                                              outfile = dict_config_setup[step][key_figure]["outfile"], 
                                              use_same_random = True, 
                                              step = step
                                              )
                addition = " "+parameters_dict["adjust_cols"][0]+" ("+parameters_dict["adjust_cols_dtypes"][0]+")."
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition=addition)
            except:
                print("Error in correlation plotting. Skipped.")
                
        #%%
        key_figure = "F09"   
        if dict_config_setup[step][key_figure]["plot_distribution_hist"] == "1":
            if "NORMALISATION" in step:
                normalisation_factors = pd.read_csv(dict_config_setup[step][key_figure]["normalization_factors"], sep="\t", index_col=0)
                
                plot_distribution_hist(normalisation_factors, 
                                       outfile = dict_config_setup[step][key_figure]["outfile"],
                                       baseline = 0
                                       )
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition="")
            else:
                print("Error in normalisation factor plotting. Skipped.")
    
        ### Legend: log2-transformed normalization factors
        
        #%%
        key_figure = "F10" 
        if dict_config_setup[step][key_figure]["plot_abundance_vs_completeness"] == "1":
            try:    
                plot_completeness_vs_abundance_scatter(data_processed, 
                                                       metadata_df, 
                                                       outfile = dict_config_setup[step][key_figure]["outfile"], 
                                                       gridsize = int(dict_config_setup[step][key_figure]["gridsize"]), 
                                                       n_bins = int(dict_config_setup[step][key_figure]["n_bins"])
                                                       )
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition="")
            except:
                print("Error in abundance vs. completeness plotting. Skipped.")
                
        #%%
        key_figure = "F11"
        if dict_config_setup[step][key_figure]["plot_overall_completeness"] == "1":
            try:
                plot_overall_completeness(data_processed, 
                                          outfile = dict_config_setup[step][key_figure]["outfile"]
                                          )
                log_legend(captions, key_figure, dict_config_setup[step][key_figure]["outfile"], addition="")
            except:
                print("Error in overall completeness plotting. Skipped.")