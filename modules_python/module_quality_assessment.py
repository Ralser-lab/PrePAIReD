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

### Generalize colors as an input_variable

### For all modules to check how either plot changes, adjust labels for summarization
def plot_qa_box_int_and_number(data_origin, metadata, data_processed=None, sample_type="QC", sample_type_column="Sample.Type", size = (10, 10), fontsize = 15, ylabel_int = "# Precursors", y_label_box = "log2(Precursor quantity)", x_label = "Samples", step = "", outfile = "Intensity_distribution_boxplot.png"):
    
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

    indices_origin = [index for index in metadata.index if metadata.loc[index,sample_type_column]==sample_type if index in data_origin.columns]
    y_origin = [sum(~np.isnan(data_origin[column].values)) for column in indices_origin]
    ax_numbers_origin.scatter(np.arange(len(indices_origin)),y_origin)
    ax_numbers_origin.set_xticks(np.arange(len(indices_origin)))
    ax_numbers_origin.set_xticklabels([])
    ax_numbers_origin.set_ylabel(ylabel_int, fontsize=fontsize)
    ax_numbers_origin.set_axisbelow(True)
    ax_numbers_origin.grid(axis="y")
    ax_numbers_origin.tick_params(labelsize=fontsize)
    
    if data_processed is not None:
        indices_processed = [index for index in metadata.index if metadata.loc[index,sample_type_column]==sample_type if index in data_processed.columns]
        y_processed = [sum(~np.isnan(data_processed[column].values)) for column in indices_processed]
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
        y_lim = (np.min([np.min(y_origin),np.min(y_processed)])*0.95,np.max([np.max(y_origin),np.max(y_processed)])*1.05)
    else:
        y_lim = (np.min(y_origin)*0.95,np.max(y_origin)*1.05)
        
    ax_numbers_origin.set_ylim(y_lim)
    ax_numbers_origin.set_xlim(-1,len(indices_origin))
    if data_processed is not None:
        ax_numbers_processed.set_ylim(y_lim)
        ax_numbers_processed.set_xlim(-1,len(indices_processed))

    data_origin_filtered = [data_origin[column].copy().dropna().values for column in indices_origin]
    ax_boxplot_origin.boxplot(data_origin_filtered, patch_artist=True, medianprops=dict(color="w"))
    ax_boxplot_origin.set_ylabel(y_label_box, fontsize=fontsize)
    ax_boxplot_origin.set_xticklabels([])
    ax_boxplot_origin.set_xlabel(x_label, fontsize=fontsize)
    ax_boxplot_origin.tick_params(labelsize=fontsize)
    ax_boxplot_origin.set_axisbelow(True)
    ax_boxplot_origin.grid(axis="y")

    if data_processed is not None:
        data_processed_filtered = [data_processed[column].copy().dropna().values for column in indices_processed]
        ax_boxplot_processed.boxplot(data_processed_filtered, patch_artist=True, medianprops=dict(color="w"))
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

def get_CV_log(values, base=2):
    cv = np.sqrt(base**(np.std(values)**2)-1)
    return cv

# For Imputation, Normalization, Correction, Summarization
def plot_CV_violin(data_origin, metadata, data_processed=None, group_column = "Sample.Type", size=(5,5), fontsize = 15, fontsize_title = 15, calc_widths=True, step="", outfile="CV_violinplot.png"):
    
    group_dict = {group:metadata[metadata[group_column]==group].index for group in np.unique(metadata[group_column].values)}
    cvs_lognormal_origin = [[get_CV_log(data_origin.loc[index,group_dict[group]]) for index in data_origin.index] for group in group_dict.keys()]
    
    if data_processed is not None:
        cvs_lognormal_processed = [[get_CV_log(data_processed.loc[index,group_dict[group]]) for index in data_processed.index] for group in group_dict.keys()]
    
    if calc_widths:
        widths_origin = [np.max(np.histogram(cvs_lognormal_origin[ind],bins=25, range=(np.min(cvs_lognormal_origin),np.max(cvs_lognormal_origin)))[0]) for ind in range(len(cvs_lognormal_origin))]
        if data_processed is not None:
            widths_processed = [np.max(np.histogram(cvs_lognormal_processed[ind],bins=25, range=(np.min(cvs_lognormal_processed),np.max(cvs_lognormal_processed)))[0]) for ind in range(len(cvs_lognormal_processed))]
            
            scale = max(max(widths_origin),max(widths_processed))
            widths_origin = widths_origin/scale*0.8
            widths_processed = widths_processed/scale*0.8
        else:
            widths_origin = widths_origin/max(widths_origin)*0.8
    else:
        widths_origin = 0.8
        widths_processed = 0.8
    
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
    ax.set_xticklabels(group_dict.keys(), fontsize = fontsize, rotation = 45, ha = "right")
    ax.set_ylabel("Coefficient of variation", fontsize = fontsize)
    ax.tick_params(axis="y", labelsize = fontsize)
    
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

# For Imputation, Normalization, Correction
### Coloring
# Corrections: columns used -> continuous scale
# Normalization: mean intensities and summed intensities (log) -> continuous scale
# Imputation: mean intensities and summed intensities (log) -> continuous scale
def plot_PCA(data_origin, metadata, data_processed = None, n_components = 4, size = (10,10), fontsize = 15, fontsize_title = 20, line_offset = 0.025, fontsize_label = 15, outfile = "../Out/test.png", continuous = True, column_color = None, use_metadata = False):

    ### Examples:
    ### Drift correction
    #continuous = True,
    #column_color = "Order",
    #use_metadata = False,
    
    ### Correction
    #continuous = False
    #column_color = "MS.Batch"
    #column_color = "Sample.Type"

    if continuous:
        cmap = "viridis"
    else:
        cmap = "Spectral"
    
    data_origin[data_origin==0] = np.nan
    if not use_metadata:
        column_values_origin = np.sum(~np.isnan(data_origin.values),axis=0)
    data_origin = data_origin.copy().dropna()
    #labels_origin = data_origin.columns
    data_origin_scaled = StandardScaler().fit_transform(data_origin.T)
    if use_metadata:
        column_values_origin = metadata.loc[data_origin.columns,column_color].values
    
    pca_origin = PCA()
    data_origin_transformed = pca_origin.fit_transform(data_origin_scaled)[:,:n_components]
    var_origin = pca_origin.explained_variance_ratio_[:n_components]
    
    if data_processed is not None:
        data_processed[data_processed==0] = np.nan
        if not use_metadata:
            column_values_processed = np.sum(~np.isnan(data_processed.values),axis=0)
        data_processed = data_processed.copy().dropna()
        #labels_processed = data_origin.columns
        data_processed_scaled = StandardScaler().fit_transform(data_processed.T)
        if use_metadata:
            column_values_processed = metadata.loc[data_processed.columns,column_color].values
        
        pca_processed = PCA()
        data_processed_transformed = pca_processed.fit_transform(data_processed_scaled)[:,:n_components]
        var_processed = pca_processed.explained_variance_ratio_[:n_components]    
    
    cmap_ = mpl.cm.get_cmap(cmap)
    
    if continuous:
        if data_processed is not None:
            vmin = min(min(column_values_origin),min(column_values_processed))
            vmax = max(max(column_values_origin),max(column_values_processed))
        else:
            vmin = min(column_values_origin)
            vmax = max(column_values_origin)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        values = set(column_values_origin)
        values_dict = {value:(1+2*ind)/(2*len(values)) for ind,value in enumerate(values)}
        column_values_origin = np.vectorize(values_dict.get)(column_values_origin)
        bounds = np.linspace(0,1,len(values)+1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap_.N)
        if data_processed is not None:
            column_values_processed = np.vectorize(values_dict.get)(column_values_processed)
    
    colors_origin = [cmap_(norm(el)) for el in column_values_origin]
    if data_processed is not None:
        colors_processed = [cmap_(norm(el)) for el in column_values_processed]
    
    fig = plt.figure(figsize=(n_components*2.5,n_components*2.5))
    gs = GridSpec(n_components,n_components, hspace= 0.1, wspace = 0.1)
    
    for i in range(n_components-1):
        for j in range(i+1,n_components):
            ax = fig.add_subplot(gs[j,i])  
            #ax.plot([0,i+1],[0,j+1])
            ax.scatter(data_origin_transformed[:,i],data_origin_transformed[:,j],alpha=0.5, c=colors_origin)
    
            ax.tick_params(labelsize = fontsize)
            set_aspect_equal(ax)
            if i != 0:
                ax.set_yticklabels([])
                #pass
            else:
                ax.set_ylabel("PC "+str(j+1)+" ("+"{:.1f}".format(var_origin[j]*100)+"%)", fontsize = fontsize, labelpad = 5)
                ax.yaxis.get_major_locator().set_params(integer=True) 
            if j != n_components-1:
                ax.set_xticklabels([])
                #pass
            else:
                ax.set_xlabel("PC "+str(i+1)+" ("+"{:.1f}".format(var_origin[i]*100)+"%)", fontsize = fontsize, labelpad = 5)
                ax.xaxis.get_major_locator().set_params(integer=True) 
            
            ax.set_axisbelow(True)
            ax.grid()
            
            if data_processed is not None:
                ax = fig.add_subplot(gs[i,j])
                #ax.plot([0,j+1],[0,i+1])
                ax.scatter(data_processed_transformed[:,j],data_processed_transformed[:,i],alpha=0.5, c=colors_processed)
    
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
                    #pass
                else:
                    ax.set_ylabel("PC "+str(i+1)+" ("+"{:.1f}".format(var_processed[i]*100)+"%)", fontsize = fontsize, labelpad = 5)
                    
                if i != 0:
                    ax.set_xticklabels([])
                    #pass
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
        #if column_color:
        #line = plt.Line2D((.05,.95-(.95/n_components)),
        #                  (.95,.05+(.95/n_components)+(line_offset)*(n_components-1)/n_components), color="k", linewidth=1, linestyle=":")
        line = plt.Line2D((s_min,s_max-(.9/n_components)),
                          (s_max,s_min+(.9/n_components)+(line_offset)*(n_components-1)/n_components), color="k", linewidth=1, linestyle=":")
    
        #else:
        #    line = plt.Line2D((.05,.95),(.95,.05+line_offset), color="k", linewidth=1, linestyle=":")
            
        fig.add_artist(line)
    
    cax = fig.add_subplot(gs[j,j])
    cax.set_visible(False)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap), ax=cax, aspect = 10, fraction = 0.8)
    cbar.ax.tick_params(labelsize=fontsize)
    if use_metadata:
        cbar.ax.set_ylabel(column_color, fontsize=fontsize_label, labelpad=10)
    else:
        cbar.ax.set_ylabel("# not NaN", fontsize=fontsize_label, labelpad=10)
    if not continuous:
        cbar.ax.set_yticks(list(values_dict.values()))
        cbar.ax.set_yticklabels(values_dict.keys())
    
    plt.savefig(outfile, bbox_inches="tight")
    return

