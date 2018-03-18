# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 22:32:54 2018

@author: Admin
"""

#import neo
#import quantities as pq
import numpy as np
import elephant.statistics as stats
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#import misc

# Relative path to data (chenge to where you saved them)
#path = '../data_resliced_first_part/'
#path = '../data_resliced/'
path = '../data/'
# bin size
#w= 1*pq.ms


# Select the data

# bin size of the histogram
bin_raw  = 0.001
bin_mean = 0.01
# time line of the histogram
t_raw    = 0.2
t_mean   = 4

for j in range(6):

    plt.figure()    
    data_idx = j
    
    
    # Load the data
#    block = np.load(path + 'data_resliced_with_stats_first_part{}.npy'.format(data_idx), encoding='latin1').item()
#    block = np.load(path + 'data_resliced_with_stats{}.npy'.format(data_idx), encoding='latin1').item()
    block = np.load(path + 'data{}.npy'.format(data_idx), encoding='latin1').item()
    
    
    
    # calculating for ll neurons
    
    his_raw_all=np.zeros((156,int(t_raw/bin_raw)-1))
    his_mean_all=np.zeros((156,int(t_mean/bin_mean)-1))
    for i in range (156):
    
    
        unit=block.list_units[i]
        sp=unit.spiketrains
        
        # Defining an ISI for all trials of a neuron
        ISI_trials = []
        ISI_no=0
        
        # defining matrix to store them
        
        for k in sp:
            
            # appending all trials
            ISI_trials.extend(stats.isi(k))
            
            
        # calculating the average of the ISIs and numbers of the ISIs
        his_bin_raw=np.arange(0,t_raw,bin_raw)
        his_bin_mean=np.arange(0,t_mean,bin_mean)
        ISI_mean = np.mean(ISI_trials)
        ISI_no=len(ISI_trials)
        ISI_trials_final=ISI_trials/(ISI_mean)
        # hsitogram
        
        # the histogram of the data
#        plt.subplot(3,1,1)
#        his=np.histogram(ISI_trials_final,bins=his_bin_mean)
#        y=his[0]/(ISI_no)
#        plt.plot(his[1][:-1],y)
#        plt.ylabel('densiity')
#        plt.xlabel('Time [s]')
#        plt.title('Final')
#        plt.legend(['bin=0.02'])
        
        plt.subplot(2,1,1)
        his_mean=np.histogram(ISI_trials_final,bins=his_bin_mean)
        his_mean_all[i,:]=his_mean[0]
        plt.plot(his_mean[1][:-1], his_mean[0])
        plt.ylabel('counts')
        plt.xlabel('Time [s]')
        plt.title('divided by mean')
        plt.legend(['bin'+str(bin_mean)])
        
        plt.subplot(2,1,2)
        his_raw=np.histogram(ISI_trials, bins=his_bin_raw)
        his_raw_all[i,:]=his_raw[0]
        plt.plot(his_raw[1][:-1], his_raw[0])
        plt.ylabel('counts')
        plt.xlabel('Time [s]')
        plt.title('raw ISI')
        plt.legend(['bin'+str(bin_raw)])
        # add a 'best fit' line
        plt.show()
    
    plt.suptitle('data'+str(j))
        
        
    
        
        