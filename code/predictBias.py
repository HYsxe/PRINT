import os
import sys
import h5py
import pandas as pd
import numpy as np
import tqdm
import pickle
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp
import scipy.stats as ss
from datetime import datetime
from keras.models import load_model
from keras.models import Sequential
from keras.layers import *
from keras.models import Model
from keras import backend as K

# Receive arguments from R
args = pd.read_csv("args.txt", header = None)
main_dir = args.values[0][0]
data_dir = args.values[1][0]
context_radius = int(args.values[2][0])
n_jobs = int(args.values[3][0])
chunk_size = int(args.values[4][0])
os.system("rm args.txt")

#########################
# Functions we will use #
#########################

# One-hot encoding of a DNA sequence.
# Input: 
# (1) seq: a string of length L consisting of A/C/G/T
# Returns: 
# (1) onehot: L-by-4 encoded matrix
def onehot_encode(seq):
    mapping = pd.Series(index = ["A", "C", "G", "T"], data = [0, 1, 2, 3])
    bases = [base for base in seq]
    base_inds = mapping[bases]
    onehot = np.zeros((len(bases), 4))
    onehot[np.arange(len(bases)), base_inds] = 1
    return onehot

# Running one-hot encoding along the sequence of a genomic region
# Starting from the left, encode every consecutive sub-sequence with length = 2 * context_radius + 1
# with step size of 1 until we reach the other end of the sequence.
# Input: 
# (1) region_seq, a string of the DNA sequence in the region of interest. Must be longer than 2 * context_radius + 1
# (2) context_radius, radius of every sub-sequence
# Returns:
# (1) one-hot encoded sequence contexts for each position in the region
def region_onehot_encode(region_seq, context_radius = 50):
    
    # Calculate length of every local sub-sequence
    context_len = 2 * context_radius + 1
    
    # If region width is L, then region_seq should be a string of length L + 2 * context_len
    region_width = len(region_seq) - 2 * context_radius
    
    if "N" in region_seq:
        return np.zeros((region_width, 4))
    else:
        # First encode the whole region sequence 
        # This prevents repetitive computing for overlapping sub-sequences
        region_onehot = np.array(onehot_encode(region_seq))
        
        # Retrieve encoded sub-sequences by subsetting the larger encoded matrix
        region_onehot = np.array([region_onehot[i : (i + context_len), :] for i in range(region_width)])
        
        return region_onehot
    
##############################################
# Prepare training, validation and test data #
##############################################

if not os.path.exists(main_dir + "/data/shared/Tn5_NN_model.h5"):
    
    if not os.path.isdir(main_dir + "/data/shared"):
        os.makedirs(main_dir + "/data/shared")
    
    # Load ground truth Tn5 bias data
    print("Loading ground truth data")
    bias_data = pd.read_csv(main_dir + "data/BAC/obsBias.tsv", sep = "\t")
    
    # One-hot encoding of sequence contexts
    print("One-hot encoding of sequence contexts")
    seqs = bias_data.loc[:, "context"].values
    with mp.Pool(n_jobs) as pool:
       onehot_seqs = list(tqdm.tqdm(pool.imap(onehot_encode, seqs), total=len(seqs)))
    onehot_seqs = np.array(onehot_seqs)

    # Transform target values to facilitate training
    target = bias_data.loc[:, "obsBias_all"].values
    target = np.log10(target + 0.01)
    target = (target / 2) + 0.5
    
    # Get the indices of all BAC regions
    mapped_regions = bias_data.loc[:, "BACInd"].values
    n_regions = len(np.unique(mapped_regions))

    # Divide all the BAC regions into training, validation and test
    # This way, data from each region only falls within one of these datasets
    print("Splitting all data into training/validation/test")
    np.random.seed(42)
    region_inds = np.unique(mapped_regions)
    np.random.shuffle(region_inds)
    training_region_inds = region_inds[:int(n_regions * 0.8)]
    val_region_inds = region_inds[int(n_regions * 0.8):int(n_regions * 0.9)]
    test_region_inds = region_inds[int(n_regions * 0.9):]

    # Find the sequence indices of the corresponding regions
    training_inds = np.array([i for i in range(len(mapped_regions)) if \
                     mapped_regions[i] in training_region_inds])
    val_inds = np.array([i for i in range(len(mapped_regions)) if \
                mapped_regions[i] in val_region_inds])
    test_inds = np.array([i for i in range(len(mapped_regions)) if \
                 mapped_regions[i] in test_region_inds])

    # Split the encoded sequences and target values into training, validation and test
    training_data = onehot_seqs[training_inds]
    training_target = target[training_inds]
    val_data = onehot_seqs[val_inds]
    val_target = target[val_inds]
    test_data = onehot_seqs[test_inds]
    test_target = target[test_inds]

    # Up-sample intervals with scarce data. 
    # This is to make sure the training is not heavily influenced with sequences with certain ranges of biases
    print("Up-sampling training data")
    n_bins = 5
    bin_width = (max(training_target) - min(training_target)) / n_bins
    bin_starts = np.arange(n_bins) * bin_width + min(training_target)
    bin_ends = np.arange(1, n_bins + 1) * bin_width + min(training_target)
    bin_inds = [(training_target > bin_starts[i]) & (training_target < bin_ends[i]) 
            for i in range(n_bins)]
    bin_n = [sum(bin_ind) for bin_ind in bin_inds]
    training_inds = np.concatenate([np.random.choice(training_inds[bin_inds[i]], max(bin_n), replace = True)
                                   for i in range(n_bins)])
    np.random.shuffle(training_inds)
    training_data = onehot_seqs[training_inds]
    training_target = target[training_inds]
    
    # Data augmentation by taking the reverse complement sequence
    training_reverse_complement = np.flip(np.flip(training_data, axis = 1), axis = 2)
    training_data = np.concatenate([training_data, training_reverse_complement])
    training_target = np.concatenate([training_target, training_target])
    np.shape(training_data)

    # Randomly shuffle augmented data
    inds = np.arange(np.shape(training_data)[0])
    np.random.shuffle(inds)
    training_data = training_data[inds]
    training_target = training_target[inds]

#################################
# Model definition and training #
#################################

if not os.path.exists(main_dir + "/data/shared/Tn5_NN_model.h5"):
    
    # Model definition
    print("Training Tn5 bias model")
    inputs = Input(shape = (np.shape(test_data[0])))
    conv_1 = Conv1D(32, 5, padding = 'same', activation = 'relu', strides = 1)(inputs)
    maxpool_1 = MaxPooling1D()(conv_1)
    conv_2 = Conv1D(32, 5, padding = 'same', activation = 'relu', strides = 1)(maxpool_1)
    maxpool_2 = MaxPooling1D()(conv_2)
    conv_3 = Conv1D(32, 5, padding = 'same', activation = 'relu', strides = 1)(maxpool_2)
    maxpool_3 = MaxPooling1D()(conv_3)
    flat = Flatten()(maxpool_3)
    fc = Dense(32,activation = "relu")(flat)
    out = Dense(1,activation = "linear")(fc)
    model = Model(inputs=inputs,outputs=out)  
    model.summary()
    model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mse'])

    # Model training
    mse = tf.keras.losses.MeanSquaredError()
    prev_loss = np.inf
    for n_epoch in range(3):

        # New training epoch
        model.fit(training_data, training_target, batch_size=64, epochs = 1, 
                  validation_data=(val_data, val_target))  

        # Get MSE loss on the valdation set after current epoch
        val_pred = np.transpose(model.predict(val_data))[0]
        mse_loss = mse(val_target, val_pred).numpy()
        print("MSE on validation set " + str(mse_loss))

        # Get pred-target correlation on the validation set after current epoch
        print("Pred-target correlation " + str(ss.pearsonr(val_target, val_pred)[0]))

        # If loss on validation set starts to increase, stop training and adopt the previous saved version
        if mse_loss > prev_loss:
            break
        else:
            prev_loss = mse_loss

            # Save current model version
            model.save(main_dir + "/data/shared/Tn5_NN_model.h5")
            
    # Load last Tn5 bias model, which should be the currently best model
    model = load_model(main_dir + "/data/shared/Tn5_NN_model.h5")

    # Model evaluation on the test set
    print("Evaluating performance on the test set")
    plt_ind = np.random.choice(np.arange(len(test_data)), 10000)
    test_pred = np.transpose(model.predict(test_data))[0]
    test_target_rev = np.power(10, (test_target - 0.5) * 2) - 0.01
    test_pred_rev = np.power(10, (test_pred - 0.5) * 2) - 0.01
    
    # Plot test results
    matplotlib.use('Agg')
    plt.figure(dpi = 100)
    plt.scatter(test_target[plt_ind], test_pred[plt_ind], s = 1)
    plt.xlabel("Target label")
    plt.ylabel("Target prediction")
    plt.title("Pearson correlation = " + str(ss.pearsonr(test_target, test_pred)[0]))
    print("Pearson correlation = " + str(ss.pearsonr(test_target, test_pred)[0]))
    plt.savefig(main_dir + "/data/shared/model_testing.pdf")
    
###############################################
# Use the model to predict Tn5 bias for regions #
###############################################

if not os.path.exists(data_dir + "predBias.h5"):
    
    # Load Tn5 bias model
    model = load_model(main_dir + "/data/shared/Tn5_NN_model.h5")
    
    # Load sequences of regions
    print("Loading sequences of regions")
    region_seqs = pd.read_csv(data_dir + "regionSeqs.txt", header = None)
    region_seqs = np.transpose(region_seqs.values)[0]
    
    # Specify the radius of sequence context and the width of the region
    context_len = 2 * context_radius + 1
    region_width = len(region_seqs[0]) - 2 * context_radius
    
    # To reduce memory usage, we chunk the region list in to smaller chunks
    starts = np.arange(len(region_seqs), step = chunk_size)
    if len(starts) > 1:
        starts = starts[:(len(starts) - 1)]
        ends = starts + chunk_size
        ends[len(ends) - 1] = len(region_seqs) + 1
    else:
        starts = [0]
        ends = [len(region_seqs) + 1]
    
    # Create folder for storing intermediate results for each chunk
    os.system("mkdir " + data_dir + "chunked_bias_pred")
    
    # Go through all chunks and predict Tn5 bias
    print("Predicting Tn5 bias for regions")
    for i in tqdm.tqdm(range(len(starts))):
    
        if os.path.exists(data_dir + "chunked_bias_pred/chunk_" + str(i) + ".pkl"):
            continue
    
        print("Processing chunk No." + str(i) + " " + datetime.now().strftime("%H:%M:%S"))
        chunk_seqs = region_seqs[starts[i]:ends[i]]
        
        # Encode sequences in the current chunk into matrices using one-hot encoding
        print("Encoding sequence contexts")
        with mp.Pool(2) as pool:
           chunk_onehot = list(tqdm.tqdm(pool.imap(region_onehot_encode, chunk_seqs), total = len(chunk_seqs)))
        
        # Use neural network to predict bias
        # For sequences containing "N", the encoded matrix will be a empty zero matrix,
        # in such cases we assign 1s as bias values
        print("Predicting bias")
        pred_bias = np.array([np.transpose(model.predict(chunk_onehot[j]))[0] \
                              if (np.sum(chunk_onehot[j]) > 1) else np.ones(region_width) \
                              for j in tqdm.tqdm(range(len(chunk_onehot)), total = chunk_size)])
        
        # Reverse transform the predicted values to the original scale
        pred_bias = np.power(10, (pred_bias - 0.5) * 2) - 0.01
    
        # Save intermediate result
        with open(data_dir + "chunked_bias_pred/chunk_" + str(i) + ".pkl", "wb") as f:
            pickle.dump(pred_bias, f)
    
    # Integrate results from all chunks
    pred_bias_all = None
    for i in tqdm.tqdm(range(len(starts))):
    
        # Load intermediate result
        with open(data_dir + "chunked_bias_pred/chunk_" + str(i) + ".pkl", "rb") as f:
            pred_bias = pickle.load(f)
    
        # Reverse transform the predicted values to the original scale
        if pred_bias_all is None:
            pred_bias_all = pred_bias
        else:
            pred_bias_all = np.concatenate([pred_bias_all, pred_bias])
    
    # Write results of the current batch to file
    output_path = data_dir + "predBias.h5"
    if os.path.isfile(output_path):
        os.system("rm " + output_path)
    hf = h5py.File(output_path, 'w')
    hf.create_dataset("predBias", data = pred_bias_all)
    hf.close()
    
    # Remove intermediate files
    os.system("rm -r " + data_dir + "chunked_bias_pred")
