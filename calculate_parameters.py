#%% SETUP AND SLICING OF THE DATA
# -*- coding: utf-8 -*-
"""
Calculates "state parameters" (= CV, CV2 etc.) for each unit and compares it to the reaction time

TODO: Add mean firing rate ( elephant.statistics.mean_firing_rate )
"""

#%% SETUP

import neo
#import quantities as pq
import numpy as np
import elephant as el
#import matplotlib.pyplot as plt
#import os as os
import useful_tools as ut

# Relative path to data (chenge to where you saved them)
path = '../data/'
resultpath = '../data_resliced/'

# Select the data


#%% Reslice trials to waiting time
for i in range(6):
    block = np.load(path + 'data{}.npy'.format(i), encoding='latin1').item()
    
    block_sliced = neo.Block()
    
    for idx, trial in enumerate(block.segments):  # for each trial
            
        # Find the index of our event in the event substructure
        on_num      = trial.events[0].annotations['trial_event_labels'].index(b'CUE-OFF')
        off_num     = trial.events[0].annotations['trial_event_labels'].index(b'GO-ON')
        
        # Find the associated times
        on_time     = trial.events[0].annotations['signal'][on_num]
        off_time    = trial.events[0].annotations['signal'][off_num]
        
        # Reslice the existing spiketrains, put them into a new neo.Segment
        seg_sliced  = neo.Segment()
        for idx_s, spike in enumerate(trial.spiketrains):
            seg_sliced.spiketrains.append(spike.time_slice(on_time, off_time))
    
        # ... and stuff the result into our new block
        block_sliced.segments.append(seg_sliced)
        block_sliced.segments[-1].annotations           = trial.annotations
        block_sliced.segments[-1].annotations['RT']     = trial.events[0].annotations['signal'][trial.events[0].annotations['trial_event_labels'].index(b'SR')] - trial.events[0].annotations['signal'][trial.events[0].annotations['trial_event_labels'].index(b'GO-ON')]

        del trial.events[0].annotations['signal']       # fixes a neo bug
        block_sliced.segments[-1].events                = trial.events
        
    ut.add_channel_and_units_v2(block_sliced)
        
    
    #%% Add state parameters to every spiketrain
    cv_per_trial    = [];
    cv2_per_trial   = [];
    for idx, trial in enumerate(block_sliced.segments):  # for each trial and spiketrain    
        for idx_st, spiketrain in enumerate(trial.spiketrains):
            if len(spiketrain) > 0:
                # CV
                trial.spiketrains[idx_st].annotations['CV']          = el.statistics.cv(spiketrain)
            
                temp_isi = el.statistics.isi(spiketrain)
                
                # CV2
                if len(temp_isi) > 1:
                    trial.spiketrains[idx_st].annotations['CV2']     = el.statistics.cv2(temp_isi)
                else:
                    trial.spiketrains[idx_st].annotations['CV2']     = np.NaN   
                
                # Firing rate
                trial.spiketrains[idx_st].annotations['FR']          = el.statistics.mean_firing_rate(spiketrain)
            else:
                trial.spiketrains[idx_st].annotations['CV']          = np.NaN
                trial.spiketrains[idx_st].annotations['CV2']         = np.NaN
                trial.spiketrains[idx_st].annotations['FR']          = 0
        
    for unit in block_sliced.list_units:
         
        temp_pg = []
        temp_sg = []
        unit.annotations['FF'] = el.statistics.fanofactor(unit.spiketrains)
        
        for train in unit.spiketrains:
            if train.annotations['belong_to_trialtype'].startswith('PG'): 
                temp_pg.append(train) 
            else: 
                temp_sg.append(train) 
    
        unit.annotations['FF_PG'] = el.statistics.fanofactor(temp_pg)
        unit.annotations['FF_SG'] = el.statistics.fanofactor(temp_pg)
        
        unit.annotations['FF_avg'] = (el.statistics.fanofactor(temp_pg ) + el.statistics.fanofactor(temp_sg)) / 2
            
    np.save(resultpath + 'data_resliced_with_stats{}.npy'.format(i), block_sliced)
        
    
    

