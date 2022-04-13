import pyabf
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
# these are a set of convenience functions for working with Homeira and Lihua's files 


# this function tries to guess the gain on the response (voltage) channel
# the logic is to compare the voltage value at time 0 and compare that to the listed RMP
# the best gain minimizes that difference
CHECK_RESPONSE_GAINS = np.array([.25, .33, .5, 1, 2, 20, 25, 50, 100, 200]) # these are the valid gains that Homeira /Lihua uses
def guess_response_gain(resp_vec, stated_rmp, offset_voltage = 0):
    # this function is a bit of logic that tries to guess the gain given an rmp value
    # try to figure out gain on response channel by comparing to RMP

    abs_diff_vec = np.abs(resp_vec[0] * CHECK_RESPONSE_GAINS + offset_voltage - stated_rmp)
    best_gain_ind = np.argmin(abs_diff_vec)
    
    rmp_abs_error = abs_diff_vec[best_gain_ind]
    best_gain = CHECK_RESPONSE_GAINS[best_gain_ind]
    return(best_gain, rmp_abs_error)

# try to figure out the gain on the stimulus channel
def get_stim_gain(stim_vec):
    min_stim = np.min(stim_vec)
    if min_stim > -1:
        stim_gain = 1000
    else:
        stim_gain = 1
    return stim_gain

# parse relevant info related to stimulus, including duration, and amplitudes
def get_stim_info(abf, stim_channel_num = 1, stim_gain = 1, stim_name = 'sweepC'):
    num_sweeps = abf.sweepCount
    stim_amps = np.zeros(num_sweeps) 
    stim_start_time = None
    stim_end_time = None
    sampling_rate = int(round(1/(abf.sweepX[2] - abf.sweepX[1]))) # manually calculate the sampling rate

    for i in range(0, num_sweeps):
        abf.setSweep(i, channel=stim_channel_num)
        sampling_rate = abf.dataRate
        if stim_name == 'sweepY':
            stim_vec = np.round(abf.sweepY * stim_gain)
        else:
            stim_vec = np.round(abf.sweepC * stim_gain)
        stim_amp = stim_vec[5000]

        stim_amps[i] = round(stim_amp)
        non_zero_inds = np.where(stim_vec == stim_amp)
        stim_duration = np.shape(non_zero_inds)[1] * 1/sampling_rate
        if stim_duration == 0:
            continue
        stim_start_ind = non_zero_inds[0][0]
        stim_end_ind = non_zero_inds[0][-1]
        
        stim_start_time = abf.sweepX[stim_start_ind]
        stim_end_time = abf.sweepX[stim_end_ind]

    ret_dict = {'stim_amp_vec' : stim_amps, 'stim_duration' : stim_duration, 
                'stim_start_time' : stim_start_time, 'stim_end_time' : stim_end_time, 'num_sweeps' : num_sweeps,
               'stim_sampling_rate' : sampling_rate}
    return(ret_dict)

# gets all relevant info about stimulus, including channel, duration, etc. and returns as dictionary
def get_stim_dict(meta_row, cell_meta_df, stim_name = 'sweepC'):
    # returns path of abf file containing stim info
    # stim channel index
    # stim gain
    # other info like num sweeps, which 
    
    row = meta_row
    f = row['cell_id'].values[0]
    fn = row['full_path'].values[0]
    recorder_name = row['recorder_name'].values[0]
    abf = pyabf.ABF(fn) # loads in the abf file

    stim_info_dict = {}    

    # figure out stim channel
    stim_chan = len(abf.channelList)-1 # this seems to be generally true
    abf.setSweep(0, channel=stim_chan)

    if stim_name == 'sweepY':
        stim_vec = abf.sweepY
    else:
        stim_vec = abf.sweepC

    stim_gain = get_stim_gain(stim_vec)

    # this infers some basic info about stim amplitudes, durations, etc.
    stim_info_dict = get_stim_info(abf, stim_chan, stim_gain = stim_gain, stim_name = stim_name)
    stim_amps = stim_info_dict['stim_amp_vec']
    sampling_rate = stim_info_dict['stim_sampling_rate']
    num_sweeps = stim_info_dict['num_sweeps']
    if np.std(stim_amps) == 0 and recorder_name == 'Homeira':
        stim_chan = 0
        abf.setSweep(0, channel=stim_chan)
        stim_vec = abf.sweepC
        stim_gain = get_stim_gain(stim_vec)
        stim_info_dict = get_stim_info(abf, stim_chan, stim_gain = stim_gain)
    elif np.std(stim_amps) == 0 and recorder_name == 'Lihua':
        # logic here is that if abf file meets these criteria, we should replace the stimulus with the one from
        # a specific abf file with available info
        abf_file_name = '14617300.abf'
        curr_row = cell_meta_df.loc[cell_meta_df['cell_id'] == abf_file_name]
        row = curr_row
        return get_stim_dict(row, cell_meta_df) # woo recursion
    stim_amps = stim_info_dict['stim_amp_vec']
    if np.std(stim_amps) == 0:
        valid_stim = False
    else:
        valid_stim = True
    ret_dict = {'stim_chan' : stim_chan, 'stim_gain' : stim_gain, 'stim_path' : fn, 'valid_stim' : valid_stim}
    
    ret_dict.update(stim_info_dict)
    return ret_dict