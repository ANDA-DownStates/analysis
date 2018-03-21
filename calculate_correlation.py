#%% SETUP AND SLICING OF THE DATA
# -*- coding: utf-8 -*-
"""
Takes the resliced datasets, extracts the parameters and then plots them 
against each other
"""

#%% SETUP

import neo
import quantities as pq
import numpy as np
import elephant as el
import matplotlib.pyplot as plt
import os as os
import useful_tools as ut
import seaborn as sns
import pandas as pd
import scipy

path        = '../data_resliced/' 
#path2 = '../data/' 

#all_blocks  = [];
FFs         = []
CVs         = []
CV2s        = []
RTs         = []
FR          = []

CVs_corr_all        = []
CV2s_corr_all       = []
FRs_corr_all        = []
CVs_corrsig_all     = []
CV2s_corrsig_all    = []
FRs_corrsig_all     = []
CVs_corr_sigonly    = []
CV2s_corr_sigonly   = []
FRs_corr_sigonly    = []

FR_means = []
CV_means = []
CV2_means = []
    
for i in range(0, 6):
    i = 2
    block_sliced = np.load(path + 'data_resliced_with_stats{}.npy'.format(i)).item()
    for seg in block_sliced.segments:
        seg.events[0].annotations.pop('signal')
    
    no_trials       = len(block_sliced.segments)
    no_units        = len(block_sliced.list_units)
    RTs             = []
    trial_ids       = []
    CVs_all         = np.zeros((no_trials, no_units))
    CV2s_all        = np.zeros((no_trials, no_units))
    FRs_all          = np.zeros((no_trials, no_units))
    
    # First collect all Trial IDs and reaction times
    for idx, segment in enumerate(block_sliced.segments):
        trial_ids.append(segment.annotations['trial_id'])
        RTs.append(segment.annotations['RT'])
        
        # Then for each trial, go through all the units and gather their parameters
        unit_id         = []
        unit_counter    = 0 # ...and idx will be my trial counter
        for channel in block_sliced.channel_indexes:
            for unit in channel.units:
                cur_spiketrain              = segment.filter(targdict=[{'channel_id': channel.annotations['channel_id']},{'unit_id': unit.annotations['unit_id']}], objects=neo.SpikeTrain)[0]
                CVs_all[idx, unit_counter]  = cur_spiketrain.annotations['CV']
                CV2s_all[idx, unit_counter] = cur_spiketrain.annotations['CV2']**2
                FRs_all[idx, unit_counter]  = cur_spiketrain.annotations['FR']
                unit_counter+=1
        
    # Calculate correlations
    CVs_corr        = np.zeros((no_units))
    CVs_corrsig     = np.zeros((no_units))
    CV2s_corr       = np.zeros((no_units))
    CV2s_corrsig    = np.zeros((no_units))
    FRs_corr        = np.zeros((no_units))
    FRs_corrsig     = np.zeros((no_units))

    RTs = np.array(RTs)         # weil besser
    for unit_idx in range(CVs_all.shape[1]):
        # CV vs RT
        nonans                      = ~np.isnan(CVs_all[:, unit_idx])
        if sum(CVs_all[:, unit_idx][nonans]) > 0:
            CVs_corr[unit_idx]      = scipy.stats.pearsonr(RTs[nonans], CVs_all[:, unit_idx][nonans])[0]
            CVs_corrsig[unit_idx]   = scipy.stats.pearsonr(RTs[nonans], CVs_all[:, unit_idx][nonans])[1]
            
        else:
            CVs_corr[unit_idx]      = np.nan        # ...or better 0 ?
            CVs_corrsig[unit_idx]   = np.nan        # ...or better 0 ?
            
        # CV2 vs RT    
        nonans                      = ~np.isnan(CV2s_all[:, unit_idx])
        if sum(CV2s_all[:, unit_idx][nonans]) > 0:
            CV2s_corr[unit_idx]     = scipy.stats.pearsonr(RTs[nonans], CV2s_all[:, unit_idx][nonans])[0]
            CV2s_corrsig[unit_idx]  = scipy.stats.pearsonr(RTs[nonans], CV2s_all[:, unit_idx][nonans])[1]

        else:
            CV2s_corr[unit_idx]     = np.nan
            CV2s_corrsig[unit_idx]  = np.nan        # ...or better 0 ?
            
        # FR vs RT     
        nonans                      = ~np.isnan(FRs_all[:, unit_idx])
        if sum(FRs_all[:, unit_idx][nonans]) > 0:
            FRs_corr[unit_idx]      = scipy.stats.pearsonr(RTs[nonans], FRs_all[:, unit_idx][nonans])[0]
            FRs_corrsig[unit_idx]   = scipy.stats.pearsonr(RTs[nonans], FRs_all[:, unit_idx][nonans])[1]
        else:
            FRs_corr[unit_idx]      = np.nan        # ...or better 0 ?
            FRs_corrsig[unit_idx]   = np.nan        # ...or better 0 ?

      
    # r    
    CVs_corr_all.append(CVs_corr)
    CV2s_corr_all.append(CV2s_corr)
    FRs_corr_all.append(FRs_corr)
    
    # p values
    CVs_corrsig_all.append(CVs_corrsig) 
    CV2s_corrsig_all.append(CV2s_corrsig)
    FRs_corrsig_all.append(FRs_corrsig)
    
    # only significant rs
    CVs_corr_sigonly.append(CVs_corr[CVs_corrsig <.05])
    CV2s_corr_sigonly.append(CV2s_corr[CV2s_corrsig <.05])
    FRs_corr_sigonly.append(FRs_corr[FRs_corrsig <.05])
    
    # Gather means per neuron
    FR_means.append(np.nanmean(FRs_all, axis=0))
    CV_means.append(np.nanmean(CVs_all, axis=0))
    CV2_means.append(np.nanmean(CV2s_all, axis=0))

# Export the means for SPSS
all_means = pd.DataFrame()
for ds_idx in range(6):
    all_means = all_means.append(pd.DataFrame({ 'DataSet' : 'DS ' + str(ds_idx),
                    'Firing Rate' : FR_means[ds_idx],
                    'CV' : CV_means[ds_idx],
                    'CV2' : CV2_means[ds_idx]
                    }))
all_means.to_csv('all_means.csv')
    
# Put everything in Panda dataframes + calculate averages per dataset
all_datasets                = pd.DataFrame()
all_datasets_sigonly        = pd.DataFrame()
all_datasets_sigonly_FR     = pd.DataFrame()
all_datasets_sigonly_CV     = pd.DataFrame()
all_datasets_sigonly_CV2    = pd.DataFrame()


for ds_idx in range(len(FRs_corr_all)):
    all_datasets = all_datasets.append(pd.DataFrame({ 'DataSet' : 'DS ' + str(ds_idx),
                           'Firing Rate' : FRs_corr_all[ds_idx],
                           'CV' : CVs_corr_all[ds_idx],
                           'CV2²' : CV2s_corr_all[ds_idx]
                         }))

for ds_idx in range(len(FRs_corr_sigonly)):
    all_datasets_sigonly_FR = all_datasets_sigonly_FR.append(pd.DataFrame({ 'DataSet' : 'DS ' + str(ds_idx),
                        'Firing Rate' : FRs_corr_sigonly[ds_idx]}))
    all_datasets_sigonly_CV = all_datasets_sigonly_CV.append(pd.DataFrame({ 'DataSet' : 'DS ' + str(ds_idx),
                        'CV' : CVs_corr_sigonly[ds_idx]}))
    all_datasets_sigonly_CV2 = all_datasets_sigonly_CV2.append(pd.DataFrame({ 'DataSet' : 'DS ' + str(ds_idx),
                        'CV2' : CV2s_corr_sigonly[ds_idx]}))    
    
#%% Plot  histogram of all neurons
f, axes = plt.subplots(2,3,figsize=(15,10))

for i in range(0, 6):
    for parameter in ['Firing Rate', 'CV', 'CV2²']:
        ax = axes[i % 2, i // 2]
        sns.distplot(all_datasets.loc[all_datasets['DataSet'] == 'DS '+ str(i), [parameter]].dropna(axis = 0), ax = ax, hist = False, label = parameter)
        ax.set_ylim(0,5)
        ax.set_xlim(-1,1)
    
        ax.set_ylabel('Maximum likelihood gaussian distribution fit')
        ax.set_xlabel('Correlation coefficient')    
    
f.savefig('./Plot/All correlations FR, CV, CV2')

# Calculate mean correlation over units per dataset
mean = []
for i in range(0, 6):
    mean.append(all_datasets.loc[all_datasets['DataSet'] == 'DS '+ str(i), 'Firing Rate'].dropna(axis = 0).mean(axis = 0))
    
    
#%% Plot histogram only of signficant ones (NOT very helpful)
f, axes = plt.subplots(2,3,figsize=(15,10))

for i in range(0, 6):
        parameter = 'Firing Rate'
        ax = axes[i % 2, i // 2]
        sns.distplot(all_datasets_FR.loc[all_datasets_sigonly_FR['DataSet'] == 'DS '+ str(i), parameter].dropna(axis = 0), ax = ax, bins = 4, hist = False, rug = True, label = parameter)
        ax.set_ylim(0,6)
        ax.set_xlim(-1,1)
    
        ax.set_ylabel('Maximum likelihood gaussian distribution fit')
        ax.set_xlabel('Correlation coefficient')    

        parameter = 'CV'
        sns.distplot(all_datasets_sigonly_CV.loc[all_datasets_sigonly_CV['DataSet'] == 'DS '+ str(i), parameter].dropna(axis = 0), ax = ax, bins = 4,hist = False, label = parameter)
        ax.set_ylim(0,6)
        ax.set_xlim(-1,1)
    
        ax.set_ylabel('Maximum likelihood gaussian distribution fit')
        ax.set_xlabel('Correlation coefficient')  

        parameter = 'CV2'
        sns.distplot(all_datasets_sigonly_CV2.loc[all_datasets_sigonly_CV2['DataSet'] == 'DS '+ str(i), parameter].dropna(axis = 0), ax = ax, bins = 4,hist = False, label = parameter)
        ax.set_ylim(0,6)
        ax.set_xlim(-1,1)
    
        ax.set_ylabel('Maximum likelihood gaussian distribution fit')
        ax.set_xlabel('Correlation coefficient')  
        
f.savefig('./Plot/correlations')

#%% Calculate measures of correlation over units per dataset

# The number of significantly correlated neurons is highest in DS 2 and 5
FR_sigcount = []
for i in range(0, 6):
    FR_sigcount.append(len(all_datasets_sigonly_FR.loc[all_datasets_sigonly_FR['DataSet'] == 'DS '+ str(i),'Firing Rate']))

      
# The summed correlation coefficient of all significant correlations is highest (in the negative) in DS 2 and 5
FR_sigsum = []
for i in range(0, 6):
    FR_sigsum.append(np.sum(all_datasets_sigonly_FR.loc[all_datasets_sigonly_FR['DataSet'] == 'DS '+ str(i),'Firing Rate']))
      
FR_sigmean = []
for i in range(0, 6):
    FR_sigmean.append(np.mean(all_datasets_sigonly_FR.loc[all_datasets_sigonly_FR['DataSet'] == 'DS '+ str(i),'Firing Rate']))
         
    
CV_sigcount = []
for i in range(0, 6):
    CV_sigcount.append(len(all_datasets_sigonly_CV.loc[all_datasets_sigonly_CV['DataSet'] == 'DS '+ str(i),'CV']))
      
CV_sigsum = []
for i in range(0, 6):
    CV_sigsum.append(np.sum(all_datasets_sigonly_CV.loc[all_datasets_sigonly_CV['DataSet'] == 'DS '+ str(i),'CV']))
 
CV2_sigcount = []
for i in range(0, 6):
    CV2_sigcount.append(len(all_datasets_sigonly_CV2.loc[all_datasets_sigonly_CV2['DataSet'] == 'DS '+ str(i),'CV2']))
      
CV2_sigsum = []
for i in range(0, 6):
    CV2_sigsum.append(np.sum(all_datasets_sigonly_CV2.loc[all_datasets_sigonly_CV2['DataSet'] == 'DS '+ str(i),'CV2']))
          

f, (ax1, ax2) = plt.subplots(1,2, figsize=(20,5))

ax1.bar([0,1,2,3,4,5], FR_sigcount)
ax1.set_xticks( [0,1,2,3,4,5])
ax1.set_xticklabels(['DS 0', 'DS 1', 'DS 2', 'DS 3', 'DS 4', 'DS 5'])
ax1.get_children()[2].set_color([150/255, 0, 0]) 
ax1.get_children()[5].set_color([150/255, 0, 0]) 
ax1.set_title('FR vs. RT: Significantly correlated neurons', fontdict={'fontsize': 20})
ax1.set_xlabel('Datasets')
ax1.set_ylabel('Number of sig. correlated neurons')
ax1.set_ylim(0, 50)

ax2.bar([0,1,2,3,4,5], FR_sigsum)
ax2.set_xticks( [0,1,2,3,4,5])
ax2.set_xticklabels(['DS 0', 'DS 1', 'DS 2', 'DS 3', 'DS 4', 'DS 5'])
ax2.get_children()[2].set_color([150/255, 0, 0]) 
ax2.get_children()[5].set_color([150/255, 0, 0]) 
ax2.set_title('FR vs. RT: Direction of correlations', fontdict={'fontsize': 20})
ax2.set_xlabel('Datasets')
ax2.set_ylabel('Summed up correlation coefficients')
ax2.set_ylim(-5,5)

f.savefig('./Plot/Correlated Neurons - FR.png')



f, (ax1, ax2) = plt.subplots(1,2, figsize=(20,5))

ax1.bar([0,1,2,3,4,5], CV_sigcount)
ax1.set_xticks( [0,1,2,3,4,5])
ax1.set_xticklabels(['DS 0', 'DS 1', 'DS 2', 'DS 3', 'DS 4', 'DS 5'])
ax1.get_children()[2].set_color([150/255, 0, 0]) 
ax1.get_children()[5].set_color([150/255, 0, 0]) 
ax1.set_title('CV vs. RT: Significantly correlated neurons', fontdict={'fontsize': 20})
ax1.set_xlabel('Datasets')
ax1.set_ylabel('Number of sig. correlated neurons')
ax1.set_ylim(0, 50)

ax2.bar([0,1,2,3,4,5], CV_sigsum)
ax2.set_xticks( [0,1,2,3,4,5])
ax2.set_xticklabels(['DS 0', 'DS 1', 'DS 2', 'DS 3', 'DS 4', 'DS 5'])
ax2.get_children()[2].set_color([150/255, 0, 0]) 
ax2.get_children()[5].set_color([150/255, 0, 0]) 
ax2.set_title('CV vs. RT: Direction of correlations', fontdict={'fontsize': 20})
ax2.set_xlabel('Datasets')
ax2.set_ylabel('Summed up correlation coefficients')
ax2.set_ylim(-5,5)

f.savefig('./Plot/Correlated Neurons - CV')



f, (ax1, ax2) = plt.subplots(1,2, figsize=(20,5))

ax1.bar([0,1,2,3,4,5], CV2_sigcount)
ax1.set_xticks( [0,1,2,3,4,5])
ax1.set_xticklabels(['DS 0', 'DS 1', 'DS 2', 'DS 3', 'DS 4', 'DS 5'])
ax1.get_children()[2].set_color([150/255, 0, 0]) 
ax1.get_children()[5].set_color([150/255, 0, 0]) 
ax1.set_title('CV2 vs. RT: Significantly correlated neurons', fontdict={'fontsize': 20})
ax1.set_xlabel('Datasets')
ax1.set_ylabel('Number of sig. correlated neurons')
ax1.set_ylim(0, 50)

ax2.bar([0,1,2,3,4,5], CV2_sigsum)
ax2.set_xticks( [0,1,2,3,4,5])
ax2.set_xticklabels(['DS 0', 'DS 1', 'DS 2', 'DS 3', 'DS 4', 'DS 5'])
ax2.get_children()[2].set_color([150/255, 0, 0]) 
ax2.get_children()[5].set_color([150/255, 0, 0]) 
ax2.set_title('CV2 vs. RT: Direction of correlations', fontdict={'fontsize': 20})
ax2.set_xlabel('Datasets')
ax2.set_ylabel('Summed up correlation coefficients')
ax2.set_ylim(-5,5)

f.savefig('./Plot/Correlated Neurons - CV2')




    