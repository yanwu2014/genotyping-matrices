import sys
import os
import subprocess as sp

cluster_file = sys.argv[1]
bam_dir = sys.argv[2]

filenames = sp.check_output('ls ' + bam_dir, shell = True)
filenames = str.strip(filenames).split('\n')
filenames = [f for f in filenames if 'bam' in f]

cell_ids = []
for f in filenames:
    f = f.split('/')[-1]
    bcs = f.split('.')[0].replace('_','.')
    batch = f.split('.')[5]
    cell_ids.append(bcs + '_' + batch)

ids_dict = dict(zip(cell_ids,filenames))

print 'Filename\tCell_ID\tCluster'
with open(cluster_file) as f:
    f.readline()
    for line in f:
        line = str.strip(line)
        fields = line.split('\t')
        cell_id = fields[0]
        cluster = fields[1]
        print ids_dict[cell_id] + '\t' + cell_id + '\t' + cluster

