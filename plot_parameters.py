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

# Load a dataset to a variable called block_sliced


# Now loop over all electrodes and correlate its state parameter across all 
# trials with the reaction time

# Extract all Fano Factors
FFs = []
for unit in block_sliced.list_units:
    FFs.append(unit.annotations['FF'])
    
    
# Average CVs and CV2s across all trials of a unit
CVs = []
for unit in block_sliced.list_units:
    temp = []
    for train in unit.spiketrains:
        if not(np.isnan(train.annotations['CV'])):
            temp.append(train.annotations['CV'])
        
    CVs.append(np.mean(temp))
    
    
    
    # Average CVs and CV2s across all trials of a unit
CV2s = []
for unit in block_sliced.list_units:
    temp = []
    for train in unit.spiketrains:
        if not(np.isnan(train.annotations['CV2'])):
            temp.append(train.annotations['CV2'])
        
    CV2s.append(np.mean(temp))
    
    
    ## TODO: PLOT THE PARAMETERS