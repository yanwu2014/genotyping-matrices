### Script for plotting the distribution of UMI reads ###

import sys
import cPickle as cp
import pysam as ps
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import spacerCalling as sc

with open(sys.argv[1]) as f:
    cell_barcodes = cp.load(f)

umi_reads = []
umi_fracs = []
for cbc,cell in cell_barcodes.items():
    molTags = cell.molTags
    totalReads = cell.gbc_reads 
    for gbc,umis in molTags.items():
        for umi,count in umis.items():
            umi_fracs.append(float(count)/totalReads)
            umi_reads.append(count)

umi_reads = np.log10(umi_reads)
umi_fracs = np.log10(umi_fracs)

plt.hist(umi_reads, bins = 80, log = True)
plt.xlabel('log10(reads per UMI)')
plt.savefig('umi_reads.png', bbox_inches = 'tight')
plt.clf()

plt.hist(umi_fracs, bins = 80, log = True)
plt.xlabel('log10(read fraction per UMI)')
plt.savefig('umi_fracs.png', bbox_inches = 'tight')
plt.clf()
