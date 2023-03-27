import numpy as np
import pandas as pd
import tqdm
import h5py
import matplotlib.pyplot as plt
from keras.models import load_model
from keras.models import Sequential
from keras.layers import *
from keras.models import Model
from keras import backend as K
from Bio import SeqIO

#########################
# Functions we will use #
#########################

def onehot_encode(seq):
    mapping = pd.Series(index = ["A", "C", "G", "T", "N"], data = [0, 1, 2, 3, 4])
    bases = [base for base in seq]
    base_inds = mapping[bases]
    onehot = np.zeros((len(bases), 5))
    onehot[np.arange(len(bases)), base_inds] = 1
    return onehot[:, :4]

def peak_onehot_encode(peak_seq, context_radius = 50):
    
    context_len = 2 * context_radius + 1
    
    # If peak width is L, then peak_seq should be a string of length L + 2 * context_len
    peak_width = len(peak_seq) - 2 * context_radius
    peak_onehot = np.array(onehot_encode(peak_seq))
    peak_onehot = np.array([peak_onehot[i : (i + context_len), :] for i in range(peak_width)])
    return peak_onehot

###########################################
# Load chromosome sequences and Tn5 model #
###########################################

# Get the list of chromosome names
genome = "sacCer3"

print("Getting pre-computed Tn5 bias for " + genome)

# Load Tn5 bias model
model = load_model("../../data/shared/Tn5_NN_model.h5")

# Get genomic sequences for each chromosome
input_file = "../../data/shared/refGenomes/" + genome +".fa"
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
seq_dict = {}
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    seq_dict[name] = sequence.upper()
    
###################################
# Genome-wide Tn5 bias prediction #
###################################
    
# We go through each chromosome and predict Tn5 bias
for chrom_name in seq_dict.keys():
    
    print("Processing " + chrom_name)
    
    # Retrieve sequence for the current chromosome
    chrom_seq = seq_dict[chrom_name]
    chrom_len = len(chrom_seq)
    context_radius = 50
    chunk_size = 100000
    chrom_seq_bias = []

    # To reduce memory usage, we chunk the chromosome sequence
    starts = np.arange(chrom_len, step = chunk_size)
    if len(starts) > 1:
        starts = starts[:(len(starts) - 1)]
        ends = starts + chunk_size
        ends[len(ends) - 1] = chrom_len + 1
    else:
        starts = [0]
        ends = [chrom_len + 1]
    
    # Extract sequence for each chunk
    chunk_seqs = []
    for i in range(len(starts)):
        if i == 0:
            chunk_seqs.append("".join(["N" for i in range(context_radius)]) + chrom_seq[starts[i]:(ends[i] + 50)])
        elif i == (len(starts) - 1):
            chunk_seqs.append(chrom_seq[(starts[i] - 50):ends[i]] + "".join(["N" for i in range(context_radius)]))
        else:
            chunk_seqs.append(chrom_seq[(starts[i] - 50):(ends[i] + 50)])
    
    # One-hot encoding of sequences and bias prediction
    for chunk_seq in tqdm.tqdm(chunk_seqs):

        # Encode local sequence context for each position in the current sequence chunk
        onehot_context = peak_onehot_encode(chunk_seq, context_radius = 50)

        # Get predicted bias for the current sequence chunk
        chunk_pred_bias = np.transpose(model.predict(onehot_context))[0]

        # Reverse transform the predicted values to the original scale
        chunk_pred_bias = np.power(10, (chunk_pred_bias - 0.5) * 2) - 0.01

        chrom_seq_bias.append(chunk_pred_bias)
    
    # Integrate prediction results from all chunks
    chrom_seq_bias = np.concatenate(chrom_seq_bias, axis = 0)

    # Save results
    hf = h5py.File("../../data/shared/precomputedTn5Bias/" + genome + "Tn5Bias.h5", 'a')
    hf.create_dataset(chrom_name, data = chrom_seq_bias)
    hf.close()
