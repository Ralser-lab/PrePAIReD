#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 16:00:21 2025

@author: oliverlemke
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib as mpl
from scipy.interpolate import UnivariateSpline
from sklearn.metrics import mean_absolute_error, mean_squared_error
from scipy.stats import pearsonr, spearmanr, pointbiserialr
from matplotlib.ticker import FormatStrFormatter
import statsmodels.api as sm
from scipy.optimize import curve_fit
from sklearn.preprocessing import SplineTransformer
from matplotlib.patches import Ellipse
from scipy.stats import chi2

### Generalize colors as an input_variable

### Possible extensions
## Options to add later
## Group-wise plots
## if sample list is provided samples vs. all
## Clustermap of DIANN metrics
## Clustermap with groups (before and after), multi layer (QC vs. biological)
## PCA with explained confounders (MPCA?)

def get_CV_log(values, base=2):
    cv = np.sqrt(float(base**(np.std(values, ddof=1)**2)-1))
    return cv

def get_indices_by_type(data_origin, metadata, data_processed=None, sample_type_column=None, sample_type = None):

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

    cmap_ = mpl.cm.get_cmap(cmap)
    
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

def get_property(data, indices, metadata=None, column_metadata=None, plot_by="Correlation", data_type="continuous", shuffle=False, seed=42, take_abs = False):
    if plot_by == "Missingness":
        property_peptide = data[indices].isna().T.sum()
    elif plot_by == "Mean":
        property_peptide = data[indices].T.mean()
    elif plot_by == "Median":
        property_peptide = data[indices].T.median()
    elif plot_by == "CV":
        property_peptide = pd.Series(index=data.index, data=[get_CV_log(data.loc[index,indices]) for index in data.index])
    elif plot_by == "norm_MAE":
        property_peptide = pd.Series(index=data.index, data = [mean_absolute_error(data.loc[index,indices].fillna(data.loc[index,indices].mean()),data.loc[index,indices].mean()*np.ones(len(indices)))/data.loc[index,indices].mean() for index in data.index])
    elif plot_by == "norm_MSE":
        property_peptide = pd.Series(index=data.index, data = [mean_squared_error(data.loc[index,indices].fillna(data.loc[index,indices].mean()),data.loc[index,indices].mean()*np.ones(len(indices)))/data.loc[index,indices].mean() for index in data.index])
    elif plot_by in ["Correlation","Corr_pvalue"]:
        if metadata is None or not column_metadata:
            raise ValueError("Metadata or Metadata column or not defined.")
        metadata_values_ = [metadata.loc[data.loc[index,indices].dropna().index,column_metadata].values for index in data.index]
        if shuffle:
            np.random.seed(seed)
            for item in metadata_values_:
                np.random.shuffle(item)
        if data_type == "continuous":
            stats = [pearsonr(metadata_values_[ind],data.loc[index,data.loc[index,indices].dropna().index]) for ind,index in enumerate(data.index)]
        elif data_type == "ordinal":
            stats = [spearmanr(metadata_values_[ind],data.loc[index,data.loc[index,indices].dropna().index]) for ind,index in enumerate(data.index)]
        elif data_type == "nominal":
            stats = [pointbiserialr(metadata_values_[ind],data.loc[index,data.loc[index,indices].dropna().index]) for ind,index in enumerate(data.index)]
        else:
            raise ValueError("data_type not defined.")
        if plot_by == "Correlation":
            if take_abs:
                property_peptide = pd.Series(index=data.index, data=[abs(el[0]) for el in stats])
            else:
                property_peptide = pd.Series(index=data.index, data=[el[0] for el in stats])
        elif plot_by == "Corr_pvalue":
            property_peptide = pd.Series(index=data.index, data=[-1*np.log10(el[1]) for el in stats])
    return property_peptide

### For all modules to check how either plot changes, adjust labels for summarization
def plot_qa_box_int_and_number(data_origin, metadata, data_processed=None, sample_type=None, sample_type_column=None, size = (10, 10), fontsize = 15, ylabel_int = "# Precursors", y_label_box = "log2(Precursor quantity)", x_label = "Samples", step = "", outfile = "Intensity_distribution_boxplot.png", column_color = None, use_metadata_coloring = False, cmap=None, continuous=False):
    
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
        y_lim = (np.min([data_processed[indices_processed].min().min(),data_origin[indices_origin].min().min()])*0.95,
                 np.max([data_processed[indices_processed].max().max(),data_origin[indices_origin].max().max()])*1.05)
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
        
    if use_mask_values:
        mask.fillna("0",inplace=True)
        mask_values = list(set(pd.melt(mask)["value"].values)^set("0"))
        mask_values.sort()
    else:
        description_dict = {mask_description.loc[index,"value"]:mask_description.loc[index,"filter"] for index in mask_description.index}
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
def plot_CV_violin(data_origin, metadata, data_processed=None, group_columns = ["Sample.Type"], size=(1,5), fontsize = 15, fontsize_title = 15, calc_widths=True, step="", outfile="CV_violinplot.png", base=2, scale_max=0.95, y_lim=None, write_to_disk=False):
    
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
                    
    #cvs_lognormal_origin = [[el for el in [get_CV_log(data_origin.loc[index,[column for column in group_dict[group] if column in data_origin.columns]], base=base) for index in data_origin.index] if not np.isnan(el)] for group in group_dict.keys()]
    cvs_lognormal_origin = [[el for el in [get_CV_log(data_origin.loc[index,[column for column in group_dict[group] if column in data_origin.columns]], base=base) for index in data_origin.index]] for group in group_dict.keys()]

    if write_to_disk:
        pd.DataFrame(columns=data_origin.index,index=group_dict.keys(),data=cvs_lognormal_origin).T.to_csv(outfile[:-4]+"_CVs_origin.tsv",sep="\t")

    cvs_lognormal_origin = [[el for el in item if not np.isnan(el)] for item in cvs_lognormal_origin]
    cvs_lognormal_origin = [item if item or sum(item)>0 else [np.nan]*2 for item in cvs_lognormal_origin]
    
    if data_processed is not None:
        #cvs_lognormal_processed = [[el for el in [get_CV_log(data_processed.loc[index,[column for column in group_dict[group] if column in data_processed.columns]], base=base) for index in data_processed.index] if not np.isnan(el)] for group in group_dict.keys()]
        cvs_lognormal_processed = [[el for el in [get_CV_log(data_processed.loc[index,[column for column in group_dict[group] if column in data_processed.columns]], base=base) for index in data_processed.index]] for group in group_dict.keys()]
        
        if write_to_disk:
            pd.DataFrame(columns=data_processed.index,index=group_dict.keys(),data=cvs_lognormal_processed).T.to_csv(outfile[:-4]+"_CVs_processed.tsv",sep="\t")
        
        cvs_lognormal_processed = [[el for el in item if not np.isnan(el)] for item in cvs_lognormal_processed]
        cvs_lognormal_processed = [item if item else [np.nan]*2 for item in cvs_lognormal_processed]
    
    if y_lim:
        cvs_lognormal_origin = [[el for el in item if el>=y_lim[0] and el<=y_lim[1]] for item in cvs_lognormal_origin]
        if data_processed is not None:
            cvs_lognormal_processed = [[el for el in item if el>=y_lim[0] and el<=y_lim[1]] for item in cvs_lognormal_processed]
    
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
    return

def set_aspect_equal(ax):
    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()
    aspect = abs((x_lim[1]-x_lim[0])/(y_lim[1]-y_lim[0]))
    ax.set_aspect(aspect)
    return

def get_confidence_ellipse(x, y, confidence=0.95, edgecolor='black', facecolor='none', linewidth=2, **kwargs):

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

# For Imputation, Normalization, Correction
### Coloring
# Corrections: columns used -> continuous scale
# Normalization: mean intensities and summed intensities (log) -> continuous scale
# Imputation: mean intensities and summed intensities (log) and missingness -> continuous scale
def plot_PCA(data_origin, metadata, data_processed = None, n_components = 4, size = (10,10), fontsize = 15, fontsize_title = 20, line_offset = 0.025, fontsize_label = 15, outfile = "PCA.png", continuous = True, column_color = None, use_metadata_coloring = False, coloring_property = "Missingness", base=2, cmap=None, plot_confidence_ellipse=False, confidence=0.95, set_axis_lim = False, axis_offset = 0.1):

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
    return

def plot_distribution(data_origin, metadata, data_processed=None, sample_type = None, sample_type_column = None, step="", outfile = "Distribution.png", size = (10,5), nbins = 50, smoothing_factor = 1e-2, fontsize = 15,xlabel = "log2(PQ)", cmap = "Spectral"):

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
    
    distribution_origin = [np.histogram(data_origin[column],bins=nbins,range=(min_val,max_val),density=True)[0] for column in indices_origin]
    
    xs = np.linspace(min_val+(max_val-min_val)/(2*nbins),max_val-(max_val-min_val)/(2*nbins),nbins)
    xs_fine = np.linspace(min_val+(max_val-min_val)/(2*nbins),max_val-(max_val-min_val)/(2*nbins),250)
    
    cmap_ = mpl.cm.get_cmap(cmap)
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
        distribution_processed = [np.histogram(data_processed[column],bins=nbins,range=(min_val,max_val),density=True)[0] for column in indices_processed]
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

### Lin Reg per selection -> flag if continuous = False
### Title Before/After
def plot_scatter_series(data_origin, metadata, column_metadata, data_processed = None, sample_type = None, sample_type_column = None, plot_by = "Correlation", outfile = "Scatter_series.png", data_type = "continuous", top_n = 3, size = (5,3), ylabel = "log2(Prec. Abundance)", fontsize = 15, use_metadata_coloring = False, cmap = None, continuous = False, column_color = None, plot_linreg = False, base = 2, step=""):

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
        
    # Calculate property
    property_peptide = get_property(data_origin, indices_origin, metadata, column_metadata=column_metadata, plot_by=plot_by, data_type=data_type, take_abs=True)
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
    
    for ind, el in enumerate(selection):
        
        # estimate y_lim
        if data_processed is not None:
            y_lim = (np.min([np.min(data_origin.loc[el,indices_origin]),np.min(data_processed.loc[el,indices_processed])])*0.98,np.max([np.max(data_origin.loc[el,indices_origin]),np.max(data_processed.loc[el,indices_processed])])*1.02)
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
        ax_scatter_origin.set_ylim(y_lim)
        ax_scatter_origin.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax_scatter_origin.set_ylabel(ylabel+"\n"+el, fontsize=fontsize)
        ax_scatter_origin.set_axisbelow(True)
        ax_scatter_origin.grid(axis="y")
        ax_scatter_origin.tick_params(labelsize=fontsize)
        
        if plot_linreg:
            indices_origin_nonan = [index for index in indices_origin if not np.isnan(metadata.loc[index,column_metadata]) and not np.isnan(data_origin.loc[el,index])]
            coef = np.polyfit(metadata.loc[indices_origin_nonan,column_metadata], data_origin.loc[el,indices_origin_nonan],1) #### NaN Handling
            x = [np.min(metadata.loc[indices_origin_nonan,column_metadata]),np.max(metadata.loc[indices_origin_nonan,column_metadata])]
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
                indices_processed_nonan = [index for index in indices_processed if not np.isnan(metadata.loc[index,column_metadata]) and not np.isnan(data_processed.loc[el,index])]
                coef = np.polyfit(metadata.loc[indices_processed_nonan,column_metadata], data_processed.loc[el,indices_processed_nonan],1) #### NaN Handling
                x = [np.min(metadata.loc[indices_processed_nonan,column_metadata]),np.max(metadata.loc[indices_processed_nonan,column_metadata])]
                linreg = np.poly1d(coef)
                ax_scatter_processed.plot(x,linreg(x),c="k")
    
    plt.savefig(outfile, bbox_inches="tight")
    return

def plot_decision_boundary_estimation_curve(data_origin, metadata, method="response-binned-splines", outfile="Decision_boundary_estimation.png", sample_type_column=None, sample_type=None, size=(5,5), n_knots=3, boundary=0.5, detection_proportion_cap=1, x_lim=None, x_label="Mean(log2(Abundance))", plot_all_fits=True, fontsize=15):
    
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
    
        hist_corr_processed = np.histogram(corr_peptide,bins=nbins,range=(range_min_processed,range_max_processed))
        if not use_same_random:
            hist_corr_random_processed = np.histogram(corr_peptide_random,bins=nbins,range=(range_min_processed,range_max_processed))
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
    try:
        metadata[column_metadata] = pd.to_numeric(metadata[column_metadata])
    except:
        raise TypeError("Meta data column needs to be numeric.")
    
    # sample_type selection
    indices_origin = get_indices_by_type(data_origin, metadata, sample_type_column, sample_type)
    
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
        if any(mask=="MAR") or any(mask=="MNAR") or any(mask=="ProtAbsence") or any("ProtAbsence/MAR") or any("ProtAbsence/MNAR"):
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