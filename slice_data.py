# -*- coding: utf-8 -*-
"""
Cuts out the waiting period from our data and adds the reaction time (RT)
"""

import neo
import quantities as pq
import numpy as np
import useful_tools as ut
     
path        = '../data/'

# For old data
# data_idx = 0    # select dataset
# block_old = np.load(path + 'data{}.npy'.format(data_idx), encoding='latin1').item() 

# For new data
block = np.load(path + 'data_lfp_py3.npy', encoding='latin1').item() # new data with lfp

# Get rid of the signal event annotations, they are a problem during later saving/loading
for seg in block.segments:
    seg.events[0].annotations.pop('signal')
     
block_sliced = neo.Block()

for idx, segment in enumerate(block.segments):  # for each trial
        
    # Find the index of our event in the event substructure
    on_num      = segment.events[0].annotations['trial_event_labels'].index(b'CUE-OFF')
    off_num     = segment.events[0].annotations['trial_event_labels'].index(b'GO-ON')
    
    # Find the associated times
    on_time     = segment.events[0][on_num]
    off_time    = segment.events[0][off_num]
    
    # Reslice the existing spiketrains, put them into a new neo.Segment
    seg_sliced  = neo.Segment()
    for spike in segment.spiketrains:
        seg_sliced.spiketrains.append(spike.time_slice(on_time, off_time))
        seg_sliced.spiketrains[-1].segment          = seg_sliced    # reroute the segment back-pointer
        seg_sliced.spiketrains[-1].unit             = None  # kill unit still pointing to the old one
    
    if len(block.segments[0].analogsignals) != 0:
        for asignal in segment.analogsignals:
            seg_sliced.analogsignals.append(asignal.time_slice(on_time, off_time))
            seg_sliced.analogsignals[-1].segment        = seg_sliced  # reroute the segment back-pointer
            seg_sliced.analogsignals[-1].channel_index  = None      # kill channel index still pointing to the old one

    # ... and stuff the result into our new block
    block_sliced.segments.append(seg_sliced)
    block_sliced.segments[-1].annotations           = segment.annotations
    block_sliced.segments[-1].events                = segment.events
    block_sliced.segments[-1].events[0].segment     = block_sliced.segments[-1]
    
ut.add_channel_and_units_v2(block_sliced)

np.save(path + 'data_lfp_py3_preptime{}.npy', block_sliced)
