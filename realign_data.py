# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 16:28:00 2018

@author: Jens
"""

import neo
import quantities as pq
import numpy as np
#import elephant as el
#import matplotlib.pyplot as plt
#import os as os
import useful_tools as ut
#import seaborn as sns
#import pandas as pd
#import scipy

path        = '../data/' 

block = np.load(path + 'data_lfp_py3.npy', encoding='latin1').item()

for seg in block.segments:
    seg.events[0].annotations.pop('signal')
    
block_sliced = neo.Block()


# find the maximal possible trial time window
earliest    = np.nan
latest      = np.nan

for segment in block.segments:
    
    onset_event    = segment.events[0].annotations['trial_event_labels'].index(b'SR')
    onset_time     = segment.events[0][onset_event]

    if onset_time.magnitude < earliest or np.isnan(earliest):
        earliest = onset_time.magnitude
    if onset_time.magnitude > latest or np.isnan(latest):
        latest = onset_time.magnitude   
    
# This is the maximum amount we can cut before and after:
cut = [earliest, np.array(segment.analogsignals[0].t_stop.magnitude//1000 - latest)] 

on_time = 0
off_time = 0
for segment in block.segments:
    # Reslice the existing spiketrains, put them into a new neo.Segment
    onset_event    = segment.events[0].annotations['trial_event_labels'].index(b'SR')
    onset_time     = segment.events[0][onset_event]
    
    seg_sliced  = neo.Segment()
    for spike in segment.spiketrains:
        seg_sliced.spiketrains.append(spike.time_slice(onset_time-cut[0]*pq.s, onset_time+cut[1]*pq.s))
        seg_sliced.spiketrains[-1].segment          = seg_sliced
        seg_sliced.spiketrains[-1].unit             = None
    
    for asignal in segment.analogsignals:
        seg_sliced.analogsignals.append(asignal.time_slice(onset_time-cut[0]*pq.s, onset_time+cut[1]*pq.s))
        seg_sliced.analogsignals[-1].segment        = seg_sliced
        seg_sliced.analogsignals[-1].channel_index  = None

    # ... and stuff the result into our new block
    block_sliced.segments.append(seg_sliced)
    block_sliced.segments[-1].annotations           = segment.annotations
    block_sliced.segments[-1].events                = segment.events
    block_sliced.segments[-1].events[0].segment     = block_sliced.segments[-1]
    

ut.add_channel_and_units_v2(block_sliced)
np.save(path + 'data_lfp_py3_resliced_new.npy', block_sliced)

