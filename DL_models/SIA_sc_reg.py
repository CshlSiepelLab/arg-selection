## Training script for inferring selection coefficient

from __future__ import print_function

import numpy as np
import pickle
import sys

from keras.models import Sequential
from keras.optimizers import Adam
from keras.layers import Dense
from keras.layers import LSTM
from keras import backend as K

import time
start = time.time()

nepoch = int(sys.argv[1]) 
dropout1 = float(sys.argv[2])
batch_size = int(sys.argv[3]) 
units1 = int(sys.argv[4]) 
units2 = int(sys.argv[5]) 
loss = sys.argv[6] 
alpha = float(sys.argv[7])
activation = sys.argv[8] 
lstm_activation = sys.argv[9] 
architecture = sys.argv[10] 

def read_input(path):
  with np.load(path) as data:
    SC_ls = data["SC"]         # shape: (N, )
    fea_Mtx = data["fea_Mtx"]  # shape: (N, 2*t+2, k)
    # N - no. of samples
    # t - no. of flanking trees in the feature extraction step
    # k - no. of discretized time points

    X = np.transpose(fea_Mtx, (0,2,1))
    Y = np.reshape(SC_ls, (-1, 1))
  return X, Y

swp_trainX, swp_trainY = read_input('TRAIN_swp.npz')
swp_valX, swp_valY = read_input('VAL_swp.npz')
swp_testX, swp_testY = read_input('TEST_swp.npz')

neu_trainX, neu_trainY = read_input('TRAIN_neu.npz')
neu_valX, neu_valY = read_input('VAL_neu.npz')
neu_testX, neu_testY = read_input('TEST_neu.npz')

trainX = np.vstack((neu_trainX,swp_trainX))
trainY = np.vstack((neu_trainY,swp_trainY))
valX = np.vstack((neu_valX, swp_valX))
valY = np.vstack((neu_valY, swp_valY))
testX = np.vstack((neu_testX, swp_testX))
testY = np.vstack((neu_testY, swp_testY))
mu = np.mean(trainX, axis=0)
sigma = np.std(trainX, axis=0)

trainX = (trainX - mu) / sigma
valX = (valX - mu) / sigma
testX = (testX - mu) / sigma

print(np.shape(trainX))
print(np.shape(testX))
print(np.shape(valX))
print(np.shape(trainY))
print(np.shape(testY))
print(np.shape(valY))

# coefficient of determination (R^2) for regression  (only for Keras tensors)
def r_square(y_true, y_pred):
  SS_res =  K.sum(K.square(y_true - y_pred))
  SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
  return ( 1 - SS_res/(SS_tot + K.epsilon()) )

print("it took", time.time() - start, "seconds.")

start = time.time()
print('Build model...')
model = Sequential()

model.add(LSTM(units1, dropout=dropout1, recurrent_dropout=dropout1, activation=lstm_activation, return_sequences=True, input_shape=(np.shape(trainX)[1], np.shape(trainX)[2])))
model.add(LSTM(units2, dropout=dropout1, recurrent_dropout=dropout1, activation=lstm_activation))
model.add(Dense(1, activation=activation))
opt = Adam(learning_rate=alpha)
model.compile(loss=loss, optimizer=opt, metrics=[r_square, 'mse', 'mae'])
print('Train...')
model.summary()
model.fit(trainX, trainY, batch_size=batch_size, epochs=nepoch, validation_data=(valX, valY))
print("it took", time.time() - start, "seconds.")

start = time.time()
pred = model.predict(testX)
print("it took", time.time() - start, "seconds.")

output = "R_" + architecture + str(nepoch) + "_predictions.txt"
np.savetxt(output, pred, fmt = '%f', delimiter=' ')
output = "R_" + architecture + str(nepoch) + "_truth.txt"
np.savetxt(output, testY, fmt = '%f', delimiter=' ')

model_name = "R_" + architecture + ".h5"
model.save(model_name) #using h5 extension