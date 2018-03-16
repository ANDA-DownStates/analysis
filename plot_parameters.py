#%% SETUP AND SLICING OF THE DATA
# -*- coding: utf-8 -*-
"""
Takes the resliced datasets, extracts the parameters and then plots them 
against each other
"""

#%% SETUP


import numpy as np

import matplotlib.pyplot as plt


# Load a dataset to a variable called block_sliced
path = '../data_resliced/'

# Defining the FFS, CV, CV2
Parameters={'data{}'.format(j): {} for j in range(6)}
plt.figure()
for i in range(6):
    block_sliced = np.load(path + 'data_resliced_with_stats{}.npy'.format(i), encoding='latin1').item()
    
    
    # Now loop over all electrodes and correlate its state parameter across all 
    # trials with the reaction time
    
    # Extract all Fano Factors
    fano=[]
    for unit in block_sliced.list_units:
        fano.append(unit.annotations['FF'])
    Parameters['data{}'.format(i)]['FFs']=fano             
            
        
    # Average CVs and CV2s across all trials of a unit
    CV = []
    for unit in block_sliced.list_units:
        temp = []
        for train in unit.spiketrains:
            if not(np.isnan(train.annotations['CV'])):
                temp.append(train.annotations['CV'])
            
        CV.append(np.mean(temp))
    Parameters['data{}'.format(i)]['CVs']=CV
        
        
        
        # Average CVs and CV2s across all trials of a unit
    CV2 =[]
    for unit in block_sliced.list_units:
        temp = []
        for train in unit.spiketrains:
            if not(np.isnan(train.annotations['CV2'])):
                temp.append((train.annotations['CV2'])**2)
            
        CV2.append(np.mean(temp))
    Parameters['data{}'.format(i)]['CV2s']=CV2
        
    
# PLOT THE PARAMETER

for i in range(6):
    plt.subplot(2,3,i+1)
    plt.plot(Parameters['data{}'.format(i)]['CV2s'],Parameters['data{}'.format(i)]['FFs'],'.')
    plt.ylabel('Fano Factor')
    plt.xlabel('$CV_2^2$')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend(['data'+str(i)])

    