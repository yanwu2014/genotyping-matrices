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

cell_barcodes_file = sys.argv[1]
minUMIReads = int(sys.argv[2])

with open(sys.argv[1]) as f:
    cell_barcodes = cp.load(f)

umi_reads = []
for k,c in cell_barcodes.items():
    totalReads = c.gbc_reads
    umis = set()
    for gbcDict in c.molTags.values():
        for umi,counts in gbcDict.items():
            if not umi in umis:
                umi_reads.append(counts)
            umis.add(umi)
    cell_barcodes[k].countBarcodes(minReads = minUMIReads)

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

sns.jointplot(x = read_fracs, y = umi_fracs, kind = 'hex', gridsize = 25,
              stat_func = None)

fig_name = sys.argv[1].replace('_cell_barcodes.pickle', '_umi_vs_read_frac.png')
plt.savefig(fig_name, bbox_inches = 'tight')
plt.clf()

fig_name = sys.argv[1].replace('_cell_barcodes.pickle', '_umi_read_dist.png')
sns.distplot(np.log(umi_reads))
plt.savefig(fig_name, bbox_inches = 'tight')
plt.clf()



