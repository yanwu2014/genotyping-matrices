import sys
import subprocess as sp
from collections import defaultdict
from scTHS_count_peaks import mergeBams
import os

# Input arguments
bam_dir = sys.argv[1]
clustering_results = sys.argv[2]
n_cores = sys.argv[3]

# Create clustering dictionary
clusters = defaultdict(list)
with open(clustering_results) as f:
    f.readline()
    for line in f:
        fields = str.strip(line).split('\t')
        filename = fields[0]
        cluster = fields[2]
        
        clusters[cluster].append(filename)

# Merge bams from the same cluster
if not os.path.isdir('cluster_merged_bams'):
    sp.call('mkdir cluster_merged_bams', shell = True)

merged_bams = []
for clust,bams in clusters.items():
    bams = [bam_dir + '/' + x for x in bams]
    
    outputBam = 'cluster_merged_bams/cluster_' + str(clust)
    merged_bams.append(outputBam)
    
    mergeBams(outputBam, bams, n_cores)

# Run MACS2 on merged bams
for bam in merged_bams:
    bam = bam + '.merged.bam'
    outFile = bam.split('.')[0].split('/')[-1]
    macsCmd = 'macs2 callpeak --nolambda --nomodel --keep-dup all -t ' + str(bam) + \
              ' -n ' + outFile + ' --outdir cluster_merged_bams'
    sp.call(macsCmd, shell = True)
