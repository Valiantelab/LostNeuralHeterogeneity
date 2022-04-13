from ipfx.sweep import Sweep, SweepSet
import pyabf
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from scripts.util_functions import guess_response_gain, get_stim_gain, get_stim_info, get_stim_dict

def cell_id_to_sweep_set(abf_file_name, meta_info_df):
    
    curr_file = abf_file_name

    meta_dict = meta_info_df

    #curr_file = '15o08020.abf'

    meta_row = meta_dict.loc[meta_dict['cell_id'] == curr_file]

    file_path = meta_row['full_path'].values[0]
    stim_file_path = meta_row['stim_path'].values[0]


    resp_abf = pyabf.ABF(file_path)
    stim_abf = pyabf.ABF(stim_file_path) # for some files we're using stim traces from a different file

    num_sweeps = int(meta_row['num_sweeps'].values[0])

    stim_channel_num = int(meta_row['stim_chan'].values[0])
    response_chan_num = int(meta_row['resp_chan'].values[0])
    stim_gain = meta_row['stim_gain'].values[0]
    response_gain = meta_row['resp_gain'].values[0]

    start_time = meta_row['stim_start_time'].values[0]
    end_time = meta_row['stim_end_time'].values[0]
    resp_sampling_rate = meta_row['resp_sampling_rate'].values[0]
    stim_sampling_rate = meta_row['stim_sampling_rate'].values[0]
    resp_offset = meta_row['resp_offset'].values[0]
    stim_name = meta_row['stim_name'].values[0]
    
    stim_dict = get_stim_info(stim_abf, stim_channel_num, stim_gain, stim_name)
    stim_amps = stim_dict['stim_amp_vec']

    # curr_epoch = (int(start_time*10000), int(end_time*10000))
    # print(curr_epoch)

    clamp_mode = "CurrentClamp"

    sweep_list = list()

    for i in range(0, num_sweeps):
        sweep_num = i
        resp_abf.setSweep(sweep_num, channel=response_chan_num)

        time_vec = resp_abf.sweepX
        response_vec = resp_abf.sweepY*response_gain + resp_offset

        stim_abf.setSweep(sweep_num, channel=stim_channel_num)
        if stim_name == 'sweepY':
            stim_vec = stim_abf.sweepY * stim_gain
        else:
            stim_vec = stim_abf.sweepC * stim_gain
        
        # sometimes, when we get stim from a different file, they have diff samp rates 0_o
        if stim_sampling_rate != resp_sampling_rate:
            new_stim_vec = np.zeros(len(time_vec))
            inds = np.where((time_vec > start_time) & (time_vec < end_time))
            new_stim_vec[inds] = stim_amps[i]
            stim_vec = new_stim_vec
            #stim_vec = signal.resample(stim_vec, len(time_vec))

        sweep = Sweep(t=time_vec,
                      v=response_vec,
                      i=stim_vec,
                      sampling_rate=resp_sampling_rate,
                      sweep_number=i,
                      clamp_mode=clamp_mode,
                      #epochs = curr_epoch
                      )
        sweep_list.append(sweep)
    sweep_set = SweepSet(sweep_list)
    return(sweep_set, start_time, end_time)