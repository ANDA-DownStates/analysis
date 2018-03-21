# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 11:56:04 2018

@author: Jens
"""

import neo
import quantities as pq
import numpy as np
import elephant as el
import matplotlib.pyplot as plt
#import os as os
import useful_tools as ut
#import seaborn as sns
#import pandas as pd
#import scipy

path        = '../data/' 

# b = np.load(path + 'data_lfp_py3_preptime.npy', encoding='latin1').item()
b = np.load(path + 'data_lfp_py3.npy', encoding='latin1').item()

for seg in b.segments:
    seg.events[0].annotations.pop('signal')
    
    
# for each unit, get all spike times
freqs   = [] # list of analyzed frequencies for coherence values
cohs    = [] # list of spike-field coherence results
sta     = [] # list of spike-triggered averages

#%% MERGE FOR VISUALIZATION (RESULT PER ELECTRODE)
#%% --------------------------------------------------------------------------
# For each channel, concatenate LFPs across trials
merged_asignals = {}
for channel in b.channel_indexes:
    merged_data = [] # channel.analogsignals[0].duplicate_with_new_array([])
    for asignal in channel.analogsignals:
        merged_data = np.append(merged_data, asignal.magnitude)
    merged_asignals[channel.analogsignals[0].duplicate_with_new_array(merged_data).annotations['connector_aligned_id']] = channel.analogsignals[0].duplicate_with_new_array(merged_data)  
    
# For each unit, concatenate spiketrains across trials    
merged_trains = []
for unit in b.list_units:
    merged_data     = []
    summed_stops    = 0 
    
    for train in unit.spiketrains:
        # We gotta delete all spikes not in our prep interval and correct the spike times by the summed time of all prior data and 
        on_time     = train.segment.events[0][train.segment.events[0].annotations['trial_event_labels'].index(b'CUE-OFF')]
        off_time    = train.segment.events[0][train.segment.events[0].annotations['trial_event_labels'].index(b'GO-ON')]
        
        temp_train  = train.time_slice(on_time, off_time).magnitude + summed_stops
        merged_data = np.append(merged_data, temp_train)    
    
        summed_stops += train.t_stop.magnitude

    temp_train = unit.spiketrains[0].duplicate_with_new_data(unit.spiketrains[0])
    temp_train.t_start = b.channel_indexes[0].analogsignals[0].t_start
    temp_train.t_stop = merged_asignals[list(merged_asignals.keys())[0]].t_stop     # take t_stop from a random merged lfp
    merged_trains.append(temp_train.duplicate_with_new_data(merged_data))
    
# For each spike train, do the calculations using a closeby LFP    
# TODO: Deal with more than one unit per connector
sta = {}  # keys: connector_aligned_ids, values: stas or coherence values
cohs = {}
counter = 1
for train in merged_trains:
    print('Starting with unit ' + str(counter))        
    # Find a closeby channel (depends entirely on the connector_aligned_id convention)
    # Find a closeby channel
    neighbor_id = train.annotations['connector_aligned_id']
    if neighbor_id < 90:
        neighbor_id += 10
        if len(b.filter(targdict = {'connector_aligned_id': neighbor_id}, objects =neo.SpikeTrain)) == 0:
            if neighbor_id % 10 == 9:
                neighbor_id -= 1
            else: 
                neighbor_id += 1
    else:
        neighbor_id -= 10
        if len(b.filter(targdict = {'connector_aligned_id': neighbor_id}, objects =neo.SpikeTrain)) == 0:
            if neighbor_id % 10 == 9:
                neighbor_id -= 1
            else: 
                neighbor_id += 1
        
    lfp = merged_asignals[neighbor_id]
    #lfp = seg.filter(targdict = {'connector_aligned_id': neighbor_id}, objects = neo.AnalogSignal)[0]
    
    train.t_stop = lfp.t_stop   # some correction of rounding errors
    train.t_start = lfp.t_start
    
    #cohs_temp, freqs                                    = el.sta.spike_field_coherence(lfp, train)
    #if len(cohs_temp) == 0:
    #    print('empty cohs!')
    #cohs[train.annotations['connector_aligned_id']]     = cohs_temp
    sta_temp                                            = el.sta.spike_triggered_average(lfp, train, (-100 * pq.ms, 100 * pq.ms))
    if len(sta_temp) == 0:
        print('empty sta!')
    sta[train.annotations['connector_aligned_id']]      = sta_temp

    counter += 1

np.save(path + 'cohs.npy', cohs)
np.save(path + 'freqs.npy', freqs)
np.save(path + 'sta.npy', sta)

#%% LFP-SPIKE COHERENCE FOR CORRELATION WITH RT
#%% --------------------------------------------------------------------------
# Collect all trial IDs
trialids = []
for seg in b.segments:
    trialids.append(seg.annotations['trial_id'])
    
# Calculate spike-field coherence for all neurons and average    
counter = 1
cohs_per_trial = {}
for trial in trialids:
    print('Starting trial ' + str(counter))
    trains = b.filter(targdict = {'trial_id_st': trial}, objects = neo.SpikeTrain)
    for train in trains:
        on_time     = train.segment.events[0][train.segment.events[0].annotations['trial_event_labels'].index(b'CUE-OFF')]
        off_time    = train.segment.events[0][train.segment.events[0].annotations['trial_event_labels'].index(b'GO-ON')]
        temp_train  = train.time_slice(on_time, off_time)     
        
        # Find a closeby channel
        neighbor_id = train.annotations['connector_aligned_id']
        if neighbor_id < 90:
            neighbor_id += 10
            if len(b.filter(targdict = {'connector_aligned_id': neighbor_id}, objects =neo.SpikeTrain)) == 0:
                if neighbor_id % 10 == 9:
                    neighbor_id -= 1
                else: 
                    neighbor_id += 1
        else:
            neighbor_id -= 10
            if len(b.filter(targdict = {'connector_aligned_id': neighbor_id}, objects =neo.SpikeTrain)) == 0:
                if neighbor_id % 10 == 9:
                    neighbor_id -= 1
                else: 
                    neighbor_id += 1
        
    # Search for the right LFP of that nearby channel and for that trial
    lfp = seg.filter(targdict = {'connector_aligned_id': neighbor_id, 'trial_id_st': trial}, objects = neo.AnalogSignal)[0]
    
    train.t_stop    = lfp.t_stop       # some correction of rounding errors
    train.t_start   = lfp.t_start
    
    cohs_temp, freqs            = el.sta.spike_field_coherence(lfp, train)
    cohs_per_trial['trial_id']  = cohs_temp
    
    counter += 1
        # sta[train.annotations['connector_aligned_id']] = el.sta.spike_triggered_average(lfp, train, (-10 * pq.ms, 10 * pq.ms))
      
        

#cohs, freqs, sta = np.load(path + 'cohs_freqs_sta.npy', encoding='latin1')
maxs = []
for coh in list(cohs.values()):
    maxs.append(max(coh))


# Plot spike-triggered average   - done
stas = []
for s in list(sta.values()):
    print(len(s))
    stas.append(s)
    
f = plt.figure(figsize = (30, 20))    
plt.plot(stas[0].times, sum(stas) / len(stas), linewidth=3.0)
plt.title('Spike-triggered LFP average\n', fontsize=48)
plt.xlabel('Time', fontsize=18)
plt.ylabel('a.u.', fontsize=18)
plt.xlim(-100, 100)
plt.axvline(x=0,color='k', linewidth=1.0)
f.savefig('Spike-triggered average across all units')
     

# Visualize electrode per electrode
frequency_range = [15, 28] # in Hz
freqs_plot = np.where(np.logical_and(freqs.magnitude >= frequency_range[0], freqs.magnitude <= frequency_range[1]))
coh_matrix = np.zeros((10,10))
for key in list(cohs.keys()):
    for f in freqs_plot[0]:
        coh_matrix[(key-1) // 10, (key-1) % 10] += cohs[key][f]
    coh_matrix[key // 10, key % 10] = coh_matrix[key // 10, key % 10] / len(list(cohs.keys()))
coh_matrix[coh_matrix == 0] = np.nan
    
fig, ax = plt.subplots(figsize = (30,30))
ra = ax.matshow(coh_matrix, cmap=plt.cm.coolwarm) #cmap=plt.cm.Blues)
ax.set_title('Coherence per electrode\n (frequency range: ' + str(frequency_range[0]) + ' - ' + str(frequency_range[1]) + ' Hz)', fontsize=48)
ax.invert_yaxis()
ra.set_clim(0,0.01)
fig.colorbar(ra)
ax.tick_params(axis = 'x', labelbottom = True)
ax.tick_params(axis = 'y', labelright = True)
ax.set_xticks(range(10))
ax.set_yticks(range(10))
fig.savefig('Spike-field coherence per electrode')
