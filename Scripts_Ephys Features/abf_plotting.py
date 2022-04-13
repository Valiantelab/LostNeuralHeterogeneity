from ipfx.sweep import Sweep, SweepSet
import pyabf
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_ephys_from_abf(abf_file_name, meta_dict, show_sweeps = [0, -1]):
    
    curr_file = abf_file_name

    #curr_file = '15o08020.abf'

    meta_row = meta_dict.loc[meta_dict['cell_id'] == curr_file]
    
    file_path = meta_row['full_path'].values[0]
    stim_file_path = meta_row['stim_path'].values[0]
    #fn = file_base_base_path + file_base_path + curr_file
    # 2016_02_04_0042.abf - example from cluster 4 - burst firing
    # 13d02049.abf - example from cluster 1

    resp_abf = pyabf.ABF(file_path)
    stim_abf = pyabf.ABF(stim_file_path) # for some files we're using stim traces from a different file

    num_sweeps = int(meta_row['num_sweeps'].values[0])
    stim_channel_num = int(meta_row['stim_chan'].values[0])
    response_chan_num = int(meta_row['resp_chan'].values[0])
    stim_gain = meta_row['stim_gain'].values[0]
    stim_name = meta_row['stim_name'].values[0]
    response_gain = meta_row['resp_gain'].values[0]
    
    stim_start_time = meta_row['stim_start_time'].values[0]
    stim_end_time = meta_row['stim_end_time'].values[0]
    resp_offset = meta_row['resp_offset'].values[0]
    
    #stim_end = 2

    sweep_num = 0
    sweep_plot_list = show_sweeps

    fig = plt.figure(figsize=(8, 5))
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)

    for i in sweep_plot_list:
        sweep_num = i
        resp_abf.setSweep(sweep_num, channel=response_chan_num)
        # plot the ADC (voltage recording)
        #ax1 = fig.add_subplot(211)
        #ax1.set_title("ADC (recorded waveform)")
        resp_vec = resp_abf.sweepY*response_gain + resp_offset
        ax1.plot(resp_abf.sweepX, resp_vec, alpha = .5)

        # plot the DAC (clamp current)
        #ax2 = fig.add_subplot(212, sharex=ax1)  # <-- this argument is new
        #ax2.set_title("DAC (stimulus waveform)")
        stim_abf.setSweep(sweep_num, channel=stim_channel_num)

        if stim_name == 'sweepY':
            stim_vec = stim_abf.sweepY * stim_gain
        else:
            stim_vec = stim_abf.sweepC * stim_gain
        #abf.setSweep(sweep_num, channel=1)
        ax2.plot(stim_abf.sweepX, stim_vec, alpha = .5)

    # decorate the plots
    ax1.set_ylabel(resp_abf.sweepLabelY)
    ax2.set_xlabel(resp_abf.sweepLabelX)
    ax2.set_ylabel(resp_abf.sweepLabelC)
    ax1.axes.set_xlim(0, stim_end_time + .25)  # <-- adjust axis like this
    return(fig, ax1, ax2)
    #plt.show()