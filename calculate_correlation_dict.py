#%% SETUP AND SLICING OF THE DATA
# -*- coding: utf-8 -*-
"""
Takes the resliced datasets, extracts the parameters and then plots them 
against each other
"""

#%% SETUP


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

no_trials       = len(block_sliced.segments)
no_units        = len(block_sliced.list_units)

# ----------------------------------------------------------

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

RTs             = []
trial_ids       = []

CV2s_all        = {}
FRs_all         = {}

#%% START
#%% ---------------------------------------------------------------------------

i = 2 # The real dataset

block_sliced = np.load(path + 'data_resliced_with_stats{}.npy'.format(i)).item()
for seg in block_sliced.segments:
    seg.events[0].annotations.pop('signal')

# First collect all Trial IDs and reaction times
for idx, segment in enumerate(block_sliced.segments):
    trial_ids.append(segment.annotations['trial_id'])
    RTs.append(segment.annotations['RT'])

for channel_idx in block_sliced.channel_indexes:
    # channel = channel_idx.annotations['channel_id']
    channel = channel_idx.units[0].spiketrains[0].annotations['connector_aligned_id']
    for unit in channel_idx.units:
       for trial in trial_ids:
           cur_spiketrain     = unit.filter(targdict={'trial_id_st': trial}, objects=neo.SpikeTrain)[0]
           
           # Get all CV2s
           if channel not in CV2s_all:
               CV2s_all[channel] = {} 
           if trial not in CV2s_all[channel]:
               CV2s_all[channel][trial] = []
           CV2s_all[channel][trial].append(cur_spiketrain.annotations['CV2']**2)
           
           # Get all FRs
           if channel not in FRs_all:
               FRs_all[channel] = {} 
           if trial not in FRs_all[channel]:
               FRs_all[channel][trial] = []
           FRs_all[channel][trial].append(cur_spiketrain.annotations['FR'])
           
           

              
# Calculate correlations
CV2s_corr       = {}
CV2s_corrsig    = {}
FRs_corr        = {}
FRs_corrsig     = {}

RTs = np.array(RTs)         # weil besser

for channel_id in list(FRs_all.keys()):
    unit_count = len(FRs_all[channel_id][5]) # 5 = random trial
    for unit_id in range(unit_count):
        FRs_for_current_unit = []
        for trial_id in trial_ids:
            FRs_for_current_unit.append(FRs_all[channel_id][trial_id][unit_id])

        temp = scipy.stats.pearsonr(RTs, FRs_for_current_unit)
        
        if channel_id not in FRs_corr:
            FRs_corr[channel_id]        = []
            FRs_corrsig[channel_id]     = []
            
        FRs_corr[channel_id].append(temp[0])
        FRs_corrsig[channel_id].append(temp[1])   
        
for channel_id, doofervalue in FRs_corr.items():
    FRs_corr[channel_id] = doofervalue[np.argmax(np.abs(doofervalue))]

for channel_id, doofervalue in FRs_corrsig.items():        # dont do this at home
    FRs_corrsig[channel_id] = doofervalue[np.argmin(np.abs(doofervalue))]


corr_matrix = np.zeros((10,10))
for key in list(FRs_corr.keys()):
    corr_matrix[(key-1) // 10, (key-1) % 10] += FRs_corr[key]

  
issig = np.zeros((10,10))
for key in list(FRs_corrsig.keys()):
    issig[(key-1) // 10, (key-1) % 10] += FRs_corrsig[key] < .05


fig, ax = plt.subplots(figsize = (30,30))
ra = ax.matshow(corr_matrix, cmap=plt.cm.coolwarm) #cmap=plt.cm.Blues)
ax.set_title('Correlations FR - RT per channel', fontsize=48)
ax.invert_yaxis()
ax.tick_params(axis = 'x', labelbottom = True)
ax.tick_params(axis = 'y', labelright = True)
ax.set_xticks(range(10))
ax.set_yticks(range(10))
ra.set_clim(-0.3,0.3)
fig.colorbar(ra)
plt.show()
fig.savefig('Correlations FR by channel')



    
fig, ax = plt.subplots(figsize = (30,30))
sig_matrix = np.array(corr_matrix)
sig_matrix[issig == 0] = np.nan
ra = ax.matshow(sig_matrix, cmap=plt.cm.coolwarm) #cmap=plt.cm.Blues)
ax.set_title('Correlations FR - RT per channel', fontsize=48)
ax.invert_yaxis()
ax.tick_params(axis = 'x', labelbottom = True)
ax.tick_params(axis = 'y', labelright = True)
ax.set_xticks(range(10))
ax.set_yticks(range(10))
ra.set_clim(-0.3,0.3)
fig.colorbar(ra)
plt.show()
fig.savefig('Correlations FR by channel - sig only.png')






