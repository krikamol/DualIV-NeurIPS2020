from __future__ import print_function

import warnings

from deepiv.models import Treatment, Response
import deepiv.architectures as architectures
import deepiv.densities as densities

import tensorflow as tf

from keras.layers import Input, Dense
from keras.models import Model
from keras.layers.merge import Concatenate

import numpy

import data_generator
import pickle

n = 5000
dropout_rate = min(1000./(1000. + n), 0.5)
epochs = int(1500000./float(n)) # heuristic to select number of epochs
epochs = 300
batch_size = 100
images = False

def datafunction(n, s, images=images, test=False):
    return data_generator.demand(n=n, seed=s, ypcor=0.5, use_images=images, test=test)

x, z, t, y, g_true = datafunction(n, 1)
data = {'x': x, 'z': z, 't': t, 'y': y, 'g_true': g_true}


print("Data shapes:\n\
Features:{x},\n\
Instruments:{z},\n\
Treament:{t},\n\
Response:{y}".format(**{'x':x.shape, 'z':z.shape,
                        't':t.shape, 'y':y.shape}))

path_to_store_data = "../saved_data/demand.pkl"
pickle.dump(data, path_to_store_data)





