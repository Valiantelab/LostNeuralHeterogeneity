from ipfx.feature_extractor import SpikeTrainFeatureExtractor, SpikeFeatureExtractor
from ipfx.stimulus_protocol_analysis import LongSquareAnalysis

from ipfx.sweep import Sweep, SweepSet
import pyabf
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import sys
from PIL import Image
import os

from scripts.util_functions import guess_response_gain, get_stim_gain, get_stim_info, get_stim_dict
from scripts.abf_plotting import plot_ephys_from_abf


# We configure a SpikeFeatureExtractor and SpikeTrainFeatureExtractor that the analysis object will use

#bessel_filter_khz = sampling_rate / 10000 * 2.5
bessel_filter_khz = 1

def get_lsa_results(sweep_set, start_time, end_time, subthresh_min_amp = -500):
    spike_extractor = SpikeFeatureExtractor(start=start_time, end=end_time, filter = bessel_filter_khz)
    spike_train_extractor = SpikeTrainFeatureExtractor(start=start_time, end=end_time, baseline_interval = .05)

    # Create the analysis object
    lsa = LongSquareAnalysis(spx=spike_extractor,
                             sptx=spike_train_extractor,
                             subthresh_min_amp= subthresh_min_amp 
                            )
    lsa_results = lsa.analyze(sweep_set)
    return lsa_results



def example_sweeps_to_fig(abf_file_name, lsa_results, cell_final_raw_meta_df):
    
    curr_file = abf_file_name
    meta_dict = cell_final_raw_meta_df
    meta_row = meta_dict.loc[meta_dict['cell_id'] == curr_file]
    recorder_name = meta_row['recorder_name'].values[0]
    layer_name = meta_row['layer_name'].values[0]
    cell_type = meta_row['cell_type'].values[0]
    # if not str(cell_type) and np.isnan(cell_type):
    #     cell_type = 'Unclassified'

    rheo_sweep_index = lsa_results['rheobase_sweep'].name
    #hero_sweep_index = lsa_results['hero_sweep'].name
    num_sweeps = len(lsa_results['sweeps'])
    max_sweep_index = np.argmax(lsa_results['spiking_sweeps']['avg_rate'])

    show_sweeps = [0, rheo_sweep_index, max_sweep_index]

    (fig, ax1, ax2) = plot_ephys_from_abf(abf_file_name, cell_final_raw_meta_df, show_sweeps)
    
    fig_title_string = abf_file_name + ', ' + recorder_name + ', ' + layer_name
    fig.suptitle(fig_title_string, fontsize=12)
    fig_file_name = 'figs/' + abf_file_name + 'sweeps.png'
    plt.savefig(fig_file_name, format='png')

def lsa_results_to_fi_fig(abf_file_name, lsa_results):

    spiking_df = lsa_results["spiking_sweeps"]
    fig = plt.figure()
    plt.plot(spiking_df["stim_amp"], spiking_df["avg_rate"], '.-')
    sns.despine()
    plt.xlabel("stimulus amplitude (pA)")
    plt.ylabel("firing rate (spikes/s)")

    fig_file_name = 'figs/' + abf_file_name + 'fi.png'
    #plt.show()
    fig.savefig(fig_file_name, format='png')
    
def get_baseball_card_fig(abf_file_name, lsa_results, meta_info_df):
    
    example_sweeps_to_fig(abf_file_name, lsa_results, meta_info_df)
    lsa_results_to_fi_fig(abf_file_name, lsa_results)
    
    image_files = ['figs/' + abf_file_name + 'sweeps.png', 
                                      'figs/' + abf_file_name + 'fi.png']
    images = [Image.open(x) for x in image_files]
    widths, heights = zip(*(i.size for i in images))

    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    for im in images:
      new_im.paste(im, (x_offset,0))
      x_offset += im.size[0]

    new_im.save('figs/' + abf_file_name + '_combined.png')
    
    [os.remove(f) for f in image_files]