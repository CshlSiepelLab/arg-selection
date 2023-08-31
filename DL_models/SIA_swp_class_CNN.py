import numpy as np
from tensorflow import keras  # use tf.keras
from tensorflow.keras.layers import (  # UpSampling2D, Conv2DTranspose
    Activation, BatchNormalization, Conv2D, Dense, Dropout, Flatten, Input,
    MaxPool2D, Reshape)
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam, RMSprop
from tensorflow.keras.utils import Sequence


def fea_scaling(taxa_cnt = 128, max_gen = 1e5):
    scale_mtx = np.empty((3, taxa_cnt-1, taxa_cnt-1))

    scale_mtx[0] = 1/taxa_cnt # F
    scale_mtx[1] = 1/max_gen # W
    scale_mtx[2] = 1/taxa_cnt # R
            
    return scale_mtx

def create_convnet():
  
  input_dims = [127, 127, 3] # channel last, use np.moveaxis() to reformat data
  inputs = Input(shape=input_dims)

  # Conv layer #1
  enc_f = Conv2D(128, [127, 3], data_format='channels_last') # kernel size = [height, width]
  h_layer = enc_f(inputs)
  h_layer = Activation('relu')(h_layer)
  h_layer = MaxPool2D(pool_size=(1, 2), padding='valid', data_format='channels_last')(h_layer)

  # Conv layer #2
  enc_f = Conv2D(256, [3, 3], data_format='channels_first', padding='valid')
  h_layer = enc_f(h_layer)
  h_layer = Activation('relu')(h_layer)
  h_layer = MaxPool2D(pool_size=(4, 4), padding='valid', data_format='channels_first')(h_layer)

  # remember dimension
  #[_, cflat, wflat, hflat] = h_layer.shape.as_list() # channels first

  h_layer = Flatten()(h_layer)

  # Dense layer #1
  # classification branch
  h_class = Dense(1024, use_bias=False)(h_layer)
  h_class = BatchNormalization()(h_class)
  h_class = Activation('relu')(h_class)
  h_class = Dropout(dropout_rate)(h_class)

  # Dense layer #2
  h_class = Dense(512, activation="relu")(h_class)
  h_class = Dropout(dropout_rate)(h_class)

  # classification output
  out_class = Dense(1, activation='sigmoid')(h_class)
  ### classification model done###
  convnet_model = Model(inputs=inputs, outputs=out_class)
  convnet_model.compile(optimizer='adam', loss='binary_crossentropy', metrics='binary_accuracy')
  ######

  return convnet_model

class XYseq(Sequence):
  def __init__(self, neu_file, swp_file, trnidx_neu, trnidx_swp, batch_size):
    self.offset = 10**6 # make sure this is way larger than # of samples in each category
    self.memmap_neu = np.load(neu_file).astype(np.int32)
    self.memmap_swp = np.load(swp_file).astype(np.int32)

    self.scaler = fea_scaling()
    self.batch_size = batch_size
    self.idx_map = np.concatenate((trnidx_neu, trnidx_swp+self.offset))
    self.dsize = len(self.idx_map)
    self.no_batch = int(np.floor(self.dsize / self.batch_size)) # model sees training sample at most once per epoch
    self.indices = np.arange(self.dsize)
    np.random.shuffle(self.indices)

  def __len__(self):
    return self.no_batch

  def __getitem__(self, idx):
    batch_idx = self.indices[idx*self.batch_size:(idx+1)*self.batch_size]
    data_idx = self.idx_map[batch_idx]

    neu_idx = data_idx[data_idx < self.offset]
    swp_idx = data_idx[data_idx >= self.offset] - self.offset

    batch_X = np.concatenate((self.memmap_neu[neu_idx], self.memmap_swp[swp_idx]))
    batch_X = batch_X*self.scaler
    batch_Y = np.concatenate((np.zeros(len(neu_idx)), np.ones(len(swp_idx))))

    batch_X = np.moveaxis(batch_X, 1, -1) # change format to "channel-last"

    return batch_X, batch_Y
    
  def on_epoch_end(self):
    np.random.shuffle(self.indices)