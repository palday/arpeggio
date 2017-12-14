import mne
from philistine.mne import abs_threshold, retrieve
import numpy as np


trg_id = {"congruent":7,
          "intermediate":8,
          "mismatch":9,
          "first/onset":11,
          "prime/onset":18,
          "second/onset":21,
          "target/onset":27,
          "response/1":1,
          "response/2":2,
          "response/3":3,
          "response/4":4,
          "response/5":5,
          "response/space":6}

music_trg_id = { ("music/"+k):trg_id[k] for k in trg_id}
lang_trg_id = { ("lang/"+k):trg_id[k] for k in trg_id}

# layout = mne.find_layout(raw.info)
conditions = ("music","language")#,"decoding")

raw = dict()
for c in conditions:
    raw[c] = mne.io.read_raw_brainvision("Raw_data/01_{}.vhdr".format(c),preload=True)
    raw[c].rename_channels({'LEOG':'LO1','LBEOG':'IO1','REOG':'LO2','C64':'SO1','RM':'A2'})
    raw[c].set_eeg_reference([])
    raw[c] = mne.add_reference_channels(raw[c],"A1")
    raw[c].set_eeg_reference(['A1','A2'])
    raw[c].set_channel_types({'LO1':'eog','IO1':'eog','LO2':'eog','SO1':'eog',"A1":'misc',"A2":'misc'})
    raw[c].filter(0.15,30,method='fir',phase='zero',n_jobs=2,fir_window='hamming',fir_design='firwin')


# peak-to-peak rejection criteria
peak_to_peak = dict(eeg=150e-6,eog=250e-6)
# flatline rejection criteria
flatline = dict(eeg=5e-6)
# absolute threshold rejection criteria
threshold = 75e-6

epochs = dict()
for c in conditions:
    events = mne.find_events(raw[c])
    if c == 'music':
        event_id = music_trg_id
    elif c == 'language':
        event_id = lang_trg_id
    else:
        raise ValueError("you messed up")
    epochs[c] = mne.Epochs(raw[c],events=events, event_id=event_id,
                        tmin=-0.2, tmax=1.2, # epochs go from -200 to +1200ms
                        detrend=1, # linear detrending
                        baseline=None, # No baselining
                        reject=peak_to_peak, # peak-to-peak rejections
                        flat=flatline, # flatline rejection
                        on_missing='ignore', # ignore missing events
                        preload=True)

    epochs[c].drop_bad()
    bad_epoch_mask = abs_threshold(epochs[c].pick_types(eeg=True), threshold)
    epochs[c].drop(bad_epoch_mask,reason="absolute threshold")


#epochs = mne.concatenate_epochs([epochs[c] for c in epochs])

lang_erp = epochs['language']['lang/target'].average()
music_erp = epochs["music"]['music/target'].average()

lang_minus_music = mne.combine_evoked([lang_erp,music_erp], [1, -1])


# windows of interest for statistics
# units are *milliseconds* relative to time-locking event
wins = dict(baseline=(-200,0),
            P100=(50,150),
            N200=(150,250),
            N400=(300,500))


for c in epochs:
    df = retrieve(epochs[c], wins)
    df.to_csv("{}.csv".format(c))
