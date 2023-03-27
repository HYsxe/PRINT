#!/usr/bin/python

import sys
import re
import pysam
import os
import tqdm
#from collections import Counter
#from contextlib import contextmanager

bamfile = sys.argv[1]
tagname = sys.argv[2]
outbam = sys.argv[3]

print("Reading in bam file:", bamfile)
print("Looking for tag:", tagname)
print("Saving to new bam file:", outbam)

def getTag(intags, tag):
    '''
    Checks for specific tag from CLI
        '''
    for tg in intags:
        if(tag == tg[0]):
                return(tg[1])
    return("NA")

# BAM I/O
bam = pysam.AlignmentFile(bamfile, "rb")
out = pysam.AlignmentFile(outbam, "wb", template = bam)

# Loop over bam and extract the sequence
for read in tqdm.tqdm(bam.fetch()):
    tag = getTag(read.tags, tagname) # Fetch cell barcode
    tag = tag.replace("_", ",")
    read.query_name = read.query_name + '_' + tag
    out.write(read)

out.close()
