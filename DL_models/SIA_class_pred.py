#!/usr/bin/env python3

import sys # command line args: sys.argv
import gc
#import glob
from os import path
import numpy as np
#from sklearn.metrics import confusion_matrix
import tensorflow as tf
from tensorflow import keras # use tf.keras
from tensorflow.keras.models import Model
#from tensorflow.keras.utils import plot_model, to_categorical
#from tensorflow.keras.layers import Input, Dense, Dropout, Flatten, BatchNormalization, Activation
#from tensorflow.keras.layers import LSTM, Bidirectional
#from tensorflow.keras.regularizers import l2

print(tf.config.experimental.list_physical_devices('GPU'), flush=True)

NUM_HAPS = int(sys.argv[1]) # maize - 80, teo - 98
FEA_PREF = sys.argv[2]      # prefix to feature matrices
NO_THR = int(sys.argv[3])   # total number of threads
MOD_PATH = sys.argv[4]      # path to DL model
OUT_HANDLE = sys.argv[5]

def reshape(fea_df, samp_size):
  fea_df = fea_df/samp_size
  fea_reshaped = np.transpose(fea_df, (0, 2, 1))
  return fea_reshaped

model = tf.keras.models.load_model(MOD_PATH)

# Initialize output data structures
SIA_pred_df = np.empty((0, 3))
pos_arr = np.empty(0, dtype=int)
daf_arr = np.empty(0)
flpflg_arr = np.empty(0, dtype=bool)

for thr in range(NO_THR):
  exists = path.isfile(f"{FEA_PREF}_thr{thr}.npz")
  if not exists:
    print(f"Thread {thr} file not found!")
    continue

  with np.load(f"{FEA_PREF}_thr{thr}.npz") as npzF:
    pos_arr = np.concatenate((pos_arr, npzF["POS"].astype(int)))
    daf_arr = np.concatenate((daf_arr, npzF["DAF"]))
    flpflg_arr = np.concatenate((flpflg_arr, npzF["FLP"]))

    X_df = reshape(npzF["FEA"], NUM_HAPS)
  
  SIA_pred_df = np.concatenate((SIA_pred_df, model.predict(X_df)))
  print(f"Thread {thr}: feature matrix {X_df.shape}", flush=True)

print(pos_arr.shape, daf_arr.shape, flpflg_arr.shape, SIA_pred_df.shape)

np.savez_compressed(OUT_HANDLE, POS=pos_arr, DAF=daf_arr, FLP=flpflg_arr, PRED=SIA_pred_df)
