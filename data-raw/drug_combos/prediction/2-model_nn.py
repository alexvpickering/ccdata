from __future__ import division

import numpy as np
import pandas as pd
import pylab as pyplot
from lasagne import layers
from nolearn.lasagne import NeuralNet, TrainSplit, BatchIterator
from lasagne.updates import adadelta, adam, nesterov_momentum
from lasagne.nonlinearities import very_leaky_rectify
from scipy.stats import spearmanr
from sklearn.utils import shuffle
from random import choice, sample
import cPickle as pickle
import theano
import os
import gc

def plot_results(net, y1, y2, x2): 
    train_loss = np.array([i["train_loss"] for i in net.train_history_])
    valid_loss = np.array([i["valid_loss"] for i in net.train_history_])
    pyplot.plot(train_loss, linewidth=3, label="train")
    pyplot.plot(valid_loss, linewidth=3, label="valid")
    pyplot.grid()
    pyplot.legend()
    pyplot.xlabel("epoch")
    pyplot.ylabel("loss")
    pyplot.ylim(y1, y2)
    pyplot.xlim(-50, x2)
    pyplot.show()
    
            
def mean_cor(y, yhat):
    nrows = y.shape[0]
    cors = []
    
    # get correlation between each predicted and actual sample
    for i in range(nrows):
        cors.append(spearmanr(y[i, :], yhat[i, :])[0])
        
    return np.mean(cors)

def accuracy(y, yhat):
    num_equal = np.sum(np.sign(y) == np.sign(yhat))
    num_total = y.shape[0] * y.shape[1]
    
    return(num_equal / num_total)

# accuracy as a fraction of average model
def scaled_accuracy(y, yhat):
    return (accuracy(y, yhat) / 0.789634)


def scaled_cor(y, yhat):
    cor = mean_cor(y, yhat)   
    return (mean_cor(y, yhat) / 0.57595569)


def mae(y, yhat):
    return np.mean(np.absolute(y-yhat))



class FlipBatchIterator(BatchIterator):

    def transform(self, Xb, yb):
        Xb, yb = super(FlipBatchIterator, self).transform(Xb, yb)

        # number of samples and genes
        ns = Xb.shape[0]
        ng = Xb.shape[1] // 2
        
        # for random half of genes and samples, swap d1 and d2
        sw = [choice([0,1]) for _ in range(ng)]
        
        a = [ i    +(s*ng) for i, s in zip(range(ng), sw)]
        b = [(i+ng)-(s*ng) for i, s in zip(range(ng), sw)]

        cols  = a + b
        rows  = sample(range(ns), ns // 2)
        
        ind_orig = np.ix_(rows, range(ng*2))
        ind_flip = np.ix_(rows, cols)
        
        Xb[ind_orig] = Xb[ind_flip]
        
        return Xb, yb


def float32(k):
    return np.cast['float32'](k)

class AdjustVariable(object):
    def __init__(self, name, start=0.03, stop=0.001):
        self.name = name
        self.start, self.stop = start, stop
        self.ls = None

    def __call__(self, nn, train_history):
        if self.ls is None:
            self.ls = np.linspace(self.start, self.stop, nn.max_epochs)

        epoch = train_history[-1]['epoch']
        new_value = float32(self.ls[epoch - 1])
        getattr(nn, self.name).set_value(new_value)
        
class PickleModel(object):
    def __init__(self, freq=10):
        self.freq = freq

    def __call__(self, nn, train_history):
        epoch = train_history[-1]['epoch']
        
        if epoch % self.freq == 0:
            with open('net_int.pickle', 'wb') as f:
                pickle.dump(nn, f, -1)
        
        
class EarlyStopping(object):
    def __init__(self, patience=100):
        self.patience = patience
        self.best_valid = np.inf
        self.best_valid_epoch = 0
        self.best_weights = None

    def __call__(self, nn, train_history):
        current_valid = train_history[-1]['valid_loss']
        current_epoch = train_history[-1]['epoch']
        if current_valid < self.best_valid:
            self.best_valid = current_valid
            self.best_valid_epoch = current_epoch
            self.best_weights = nn.get_all_params_values()
        elif self.best_valid_epoch + self.patience < current_epoch:
            print("Early stopping.")
            print("Best valid loss was {:.6f} at epoch {}.".format(
                self.best_valid, self.best_valid_epoch))
            nn.load_params_from(self.best_weights)
            raise StopIteration()






# load data
X = np.load('data/X.npy')
y = np.load('data/y.npy')


# shuffle data
ids = shuffle(range(y.shape[0]), random_state=0)
X   = X[ids]
y   = y[ids]


gc.collect()


#train model
net = NeuralNet(
    layers=[
        ('input',   layers.InputLayer),
        ('dropout1', layers.DropoutLayer),
        ('hidden',  layers.DenseLayer),
        ('dropout2', layers.DropoutLayer),
        ('output',  layers.DenseLayer),
        ],
    # layer parameters:
    input_shape         = (None, X.shape[1]),
    dropout1_p          = 0.85,
    dropout2_p          = 0.5,
    hidden_num_units    = 2500,
    hidden_nonlinearity = very_leaky_rectify,
    output_nonlinearity = None, 
    output_num_units    = y.shape[1],

    # optimization method:
    batch_iterator_train    = FlipBatchIterator(batch_size=70),
    train_split             = TrainSplit(eval_size=0.2),
    update                  = nesterov_momentum,
    update_learning_rate    = theano.shared(float32(0.01)),
    update_momentum         = theano.shared(float32(0.9)),
    regression              = True, 
    max_epochs              = 5000,
    verbose                 = 1,
    on_epoch_finished       = [AdjustVariable('update_learning_rate', start=0.01, stop=0.00001),
                               AdjustVariable('update_momentum', start=0.9, stop=0.999)],
    custom_scores           = [("spr", lambda y, yhat: mean_cor(y, yhat)),
                               ("acc", lambda y, yhat: accuracy(y, yhat)),
                               ("mae", lambda y, yhat: mae(y, yhat))]
    )

np.random.seed(0)
net.fit(X, y)

plot_results(net, 0, 8, 4000)
plot_results(net, 4.8, 5.3, 10000)


# With FlipBatchIterator:
# -----------------------

# | Hidden (bs)| Dropout |Update(lr*10^4, m*10)| Activation | Epochs | Train | Val | Train/Val  | spr  | acc  | mae  | 
# | ---------- | ------- | ------------------- | ---------- | ------ |------ | ----| ---------- |----- | ---- | ---  |
# | 1000 (128) | 0.80_0.3|nm (300-1.00, 9-9.99)| vlReLU     | 4000   | 0.51  | 5.0 |  0.10      |0.508 |0.681 |1.248 |
# | 1000 (128) | 0.80_0.5|nm (300-1.00, 9-9.99)| vlReLU     | 3463   | 0.91  | 5.0 |  0.18      |0.505 |0.679 |1.255 |
# | 1000 (128) | 0.85_0.5|nm (300-1.00, 9-9.99)| vlReLU     | 4000   | 1.13  | 5.1 |  0.22      |0.498 |0.676 |1.266 |
# | 2000 (128) | 0.85_0.5|nm (300-1.00, 9-9.99)| vlReLU     | 2148   | 0.76  | 5.1 |  0.15      |0.505 |0.679 |1.255 |
# | 2000 (128) | 0.85_0.5|nm (100-0.10, 9-9.99)| vlReLU     | 4000   | 0.75  | 5.1 |  0.15      |0.508 |0.680 |1.251 |
# | 8000 (70)  | 0.85_0.5|nm (100-0.10, 9-9.99)| vlReLU     | 4073   | 0.50  | 5.0 |  0.09      |0.513 |0.683 |1.248 |
# | 2500 (70)  | 0.85_0.5|nm (100-0.10, 9-9.99)| vlReLU     | 5000   | 0.56  | 5.0 |  0.11      |0.511 |0.682 |1.244 |  # !
# | 3000 (70)  | 0.85_0.5|nm (200-0.01, 9-9.99)| vlReLU     | 6735   | 0.55  | 5.0 |  0.11      |0.511 |0.682 |1.243 |



#    257 samples
#    0.2 validation
#
#    257 * 0.8 
#  = 206 train samples
#
#    206 / 128 => 128, 78,           # uneven batch sizes
#    206 / 70  =>  70, 70, 70, 66    # more even
# 
#    2 models train on half each
# => 129, 128 train samples  
#
#    129 / 70 => 70, 59
#    128 / 70 => 70, 58



# Without FlipBatchIterator:
# --------------------------


# | Hidden (bs)| Dropout | Update              | Activation | Epochs | Train | Val | Train/Val  | spr  | acc  | mae  | 
# | ---------- | ------- | ------------------- | ---------- | ------ |-------| --- | ---------- |----- | ---- | ---  |
# | 1000 (128) | 0.8     |adadelta             | vlReLU     | 1000   | 0.19  | 5.5 |  0.03      |0.469 |0.666 |1.30  |
# | 1000 (128) | 0.8_0.1 |adadelta             | vlReLU     | 1000   | 0.40  | 5.6 |  0.07      |0.467 |0.664 |1.31  |
# | 1000 (128) | 0.8_0.1 |adadelta             | vlReLU     | 500    | 0.57  | 5.6 |  0.10      |0.459 |0.661 |1.32  |
# | 1000 (128) | 0.7_0.1 |adadelta             | vlReLU     | 1000   | 0.28  | 5.5 |  0.05      |0.468 |0.665 |1.30  |
# | 500  (128) | 0.7_0.1 |adadelta             | vlReLU     | 1000   | 0.36  | 5.5 |  0.06      |0.468 |0.665 |1.30  |  # !
# | 500  (128) | 0.7_0.2 |adadelta             | vlReLU     | 1000   | 0.69  | 5.5 |  0.12      |0.462 |0.662 |1.31  |
# | 250  (128) | 0.7_0.1 |adadelta             | vlReLU     | 1000   | 0.63  | 5.6 |  0.11      |0.459 |0.661 |1.32  |
# | 100  (128) | 0.7_0.1 |adadelta             | vlReLU     | 1000   | 1.15  | 5.7 |  0.20      |0.421 |0.645 |1.35  |
# | 1000 (128) | 0.9     |adadelta             | vlReLU     | 1000   | 0.528 | 5.7 |  0.091     |0.456 |0.659 |1.332 |
# | 1000 (128) | 0.2     |adadelta             | vlReLU     | 1000   | 0.07  | 5.7 |  0.01      |0.436 |0.651 |1.35  |
# | 1000 (128) | 0.4     |adadelta             | vlReLU     | 1000   | 0.09  | 5.6 |  0.01      |0.449 |0.657 |1.33  |
# | 1000 (128) | 0.6     |adadelta             | vlReLU     | 1000   | 0.12  | 5.6 |  0.02      |0.459 |0.661 |1.32  |
# | 2000 (128) | 0.7_0.1 |adadelta             | vlReLU     | 1000   | 0.27  | 5.5 |  0.05      |0.465 |0.663 |1.31  |