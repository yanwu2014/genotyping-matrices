import sys
import os
import os.path
import numpy as np
import pandas as pd
import subprocess as sp
import pybedtools as pyb
from functools import partial
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
from itertools import islice

# Helper function that counts the number of reads in the bam file mapped to the peaks
def _countBam_(bam_file, peaksBed):
    bamBed = pyb.BedTool(bam_file)
    return bamBed.intersect(b = peaksBed).count()

# Counts total number of unique reads mapped to called peaks in each bam file.
# Generates a log file and a knee plot
def filterBackgroundCells(sample_name, peaks_file, bam_dir, n_cores, expected_cells):
    # Load bed file
    peaksBed = pyb.BedTool(peaks_file).remove_invalid().saveas().fn 
    
    # Define countBam function for a given peaks annotation
    _countBamPartial_ = partial(_countBam_, peaksBed = peaksBed)

    # Get all bam files
    bam_files = os.listdir(bam_dir)
    bam_files = [bam_dir + '/' + f for f in bam_files if f.endswith('.bam')]
    
    # Define a pool of processors and count the unique reads mapping to called peaks
    # for each bam file
    pool = mp.Pool(n_cores)
    bam_counts = pool.map_async(_countBamPartial_, bam_files).get(9999999)
    
    # "Zip" unique read counts and file names, sort by read counts
    bam_counts = zip(bam_files,bam_counts)
    bam_counts = sorted(bam_counts, key = lambda x: x[1], reverse = True)
    
    # Output read counts to a summary file
    with open(sample_name + '_tag_counts_summary.txt', 'w') as f:
        for filename,count in bam_counts:
            f.write(filename + '\t' + str(count) + '\t' + '\n')
    
    counts = [x[1] for x in bam_counts]
    counts = np.cumsum(counts)
    
    # Plot 3 times the expected number of cells for the knee plot
    xlim = expected_cells*3
    counts = counts[0:xlim]

    # Create a knee plot
    plt.plot(range(0,len(counts)),counts)
    plt.savefig(sample_name + '_knee_plot.png')
    
    return bam_counts

# Helper function that quantifies number of reads mapping to each peak in the bam file
def _countFeatures_(bam_file, peaks_file):
    peaksBed = pyb.BedTool(peaks_file).remove_invalid()
    bamBed = pyb.BedTool(bam_file)
    bamBedCounted = peaksBed.coverage(bamBed, counts = True, 
                                      output = bam_file.replace('.bam','.counted.bed'))
    
    count_field = len(bamBedCounted[0].fields) - 1
    counts = pd.Series([x[count_field] for x in bamBedCounted])
    
    return counts

# Create a raw count matrix
def generateCountMatrix(sample_name, peaks_file, n_cores, files_to_use, min_tags):
    _countFeaturesPartial_ = partial(_countFeatures_, peaks_file = peaks_file)
    
    pool = mp.Pool(n_cores)
    cell_counts = pool.map_async(_countFeaturesPartial_, files_to_use).get(9999999)
    
    cellnames = []
    for f in files_to_use:
        bcs = f.split('/')[-1].split('.')[0].replace('_','.')
        batch = f.split('/')[-1].split('.')[5]
        cellnames.append(bcs + '_' + batch)

    cell_counts = {cellnames[i]:cell_counts[i] for i in range(0,len(files_to_use))}
    count_matrix = pd.DataFrame.from_dict(cell_counts)
    
    print 'Count Matrix Dimensions: ' + str(count_matrix.shape)
    count_matrix.to_csv(sample_name + '_count_matrix.tsv', sep = '\t')

# Merge all used files into bam for export into UCSC genome browser
def mergeBams(sample_name, files_to_use, n_cores):
    # Merge all bam files corresponding to used cells for output to genome browser
    #with open(sample_name + '_used_files.bamlist','w') as f:
    #        for fi in files_to_use:
    #            f.write(os.path.abspath(fi) + '\n')
    
    output_file = sample_name + '.merged.bam'
    # Split merge job to get around the argument limit
    if files_to_use > 1000:
        k = len(files_to_use)/1000 + 1
        for i in range(0,k):
            curr_files = files_to_use[i*1000:(i+1)*1000]
            fileArgs = 'temp_merged_' + str(i) + '.bam'
            for f in curr_files:
                fileArgs = fileArgs + ' ' + f

            mergeCmd = 'samtools merge --threads ' + str(n_cores) + ' ' + fileArgs
            sp.call(mergeCmd, shell = True)
        
        finalMergeCmd = 'samtools merge --threads ' + str(n_cores) + ' ' + output_file + ' temp_merged_*.bam' 
        sp.call(finalMergeCmd, shell = True)
        sp.call('rm temp_merged_*.bam', shell = True) # Clean up
    
    else:
        bamList = sample_name + '.bamlist'
        mergeCmd = 'samtools merge --threads ' + str(n_cores) + ' -b ' + bamList + ' ' + output_file 


def main():
    # Input arguments
    sample_name = sys.argv[1] # Unique sample id that will be added to all output files
    peaks_file = sys.argv[2] # Name of bed file with Dfilter peaks
    bam_dir = sys.argv[3] # Name of the directory with bam files for each barcode combination
    n_cores = int(sys.argv[4]) # Number of cores to use
    expected_cells = int(sys.argv[5]) # Expected number of cells (for knee plot purposes only)
    min_tags = int(sys.argv[6]) # Minimum number of unique reads a cell can have
    
    # Get total unique reads mapped to peaks for each bam file
    bam_counts = filterBackgroundCells(sample_name, peaks_file, bam_dir, n_cores, expected_cells)
    
    #bam_counts = []
    #with open('hR11_cer_tag_counts_summary.txt') as f:
    #    for line in f:
    #        fields = str.strip(line).split('\t')
    #        bam_counts.append((fields[0],int(fields[1])))

    # Keep cells with counts above a certain minimum read threshold
    files_to_use = [x[0] for x in bam_counts if x[1] > min_tags]
    print 'Number of files: ' + str(len(files_to_use))

    # Generate the count matrix
    generateCountMatrix(sample_name, peaks_file, n_cores, files_to_use, min_tags)
    
    # Merge all bam files corresponding to used cells for output to genome browser
    mergeBams(sample_name, files_to_use, n_cores)

if __name__ == '__main__': main()

