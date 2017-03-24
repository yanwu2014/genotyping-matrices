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

for k in cell_barcodes:
    cell_barcodes[k].countBarcodes(minFrac = 0)

read_fracs = []
umi_fracs = []

for k,c in cell_barcodes.items():
    totalReads = np.sum(c.gbcReadCounts.values())
    totalUMIs = np.sum(c.gbcUMICounts.values())
    for gbc,reads in c.gbcReadCounts.items():
        umis = c.gbcUMICounts[gbc]
        read_fracs.append(float(reads)/totalReads)
        umi_fracs.append(float(umis)/totalUMIs)

read_fracs = np.array(read_fracs)
umi_fracs = np.array(umi_fracs)

sns.jointplot(x = read_fracs, y = umi_fracs, kind = 'hex', bins = 'log', gridsize = 25,
              stat_func = None)

fig_name = sys.argv[1].replace('_cell_barcodes.pickle', '_umi_vs_read_frac.png')
plt.savefig(fig_name, bbox_inches = 'tight')
plt.clf()

