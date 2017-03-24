### Module with functions for handling dropSeq and 10X data ###

import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import sys
import cPickle as cp
from dedupBarcodes import collapseBarcodesExact,isInt
from dedup_umi import ClusterAndReducer
import gzip
import pysam as ps
import spacerCalling as sc
from itertools import islice
from collections import defaultdict,Counter
import subprocess as sp

def main():
    pass

   
# Input: dataFrame, cellBarcodes object, barcode to gene mapping
# Output: genotyped dataframe
def callBarcodes(dataFrameFile, cellBarcodes, barcodeToGeneFile, readThresh,
                 umiThresh, edit_dist = 1, minUMIFrac = 0.0, minGBCreads = 5,
                 minGBCumis = 2):
    
    with open(barcodeToGeneFile) as f:
        barcodeToGene = cp.load(f)
        
    for cell_bc in cellBarcodes:
        cellBarcodes[cell_bc].countBarcodes(minFrac = minUMIFrac)
        cellBarcodes[cell_bc].callBarcode(barcodeToGene, readThresh, umiThresh,
                                          minGBCreads, minGBCumis)
    
    noPlasmids = 0
    gbc_dist = Counter()
    for cell in cellBarcodes.values():
        if cell.type == 'usable':
            gbc_dist[len(cell.genotype.split(','))] += 1
        elif cell.type == 'noGRNA':
            gbc_dist[0] += 1
        else:
            noPlasmids += 1

    df = pd.read_csv(dataFrameFile, sep = '\t', header = 0, index_col = 0)
    genotyped_df = callGenotype(df, cellBarcodes)
    name = dataFrameFile.replace('.tsv','.genotyped.tsv')
    genotyped_df.to_csv(name, sep = '\t')

    print 'Total_Cells\t' + str(len(cellBarcodes))
    print 'No_Plasmid\t' + str(noPlasmids)
    
    with open('moi_distribution.txt', 'w') as f:
        for k,v in gbc_dist.items():
            f.write(str(k) + '\t' + str(v) + '\n')
    
    return cellBarcodes


# Input: DropSeq dataframe output, cellBarcode output from getCellBarcodes()
# Output: dataframe with genotypes appended after cell names
def callGenotype(dataFrame, cellBarcodes):
    dataFrameCells = set(dataFrame.columns)
    cellKeys = [key for key,cell in cellBarcodes.items() if key in dataFrameCells \
                and cell.type == 'usable']
    
    dataFrame = dataFrame[cellKeys]
    newColNames = [barcode + '_' + cellBarcodes[barcode].type + '_' + \
                   cellBarcodes[barcode].genotype for barcode in dataFrame.columns]
    # sanity check
    for i,item in enumerate(newColNames):
        assert item.split('_')[0] == dataFrame.columns[i]

    dataFrame.columns = newColNames
    return dataFrame


# Input: data frame, desired edit distance to be collapsed
# Output: collapsed data frame where cell barcodes within edit_dist have their count arrays added 
def collapseCellBarcodesDataFrame(df, edit_dist, minCounts = 1000):
    genes = list(df.index)
    df = df.to_dict('series')
    dfSums = {k:np.sum(v) for k,v in df.items()}
    cs = ClusterAndReducer()

    print 'Original_cells\t' + str(len(df))
    
    adj_list = cs.get_adj_list(df.keys(), dfSums, edit_dist, countRatio = 1.25)
    clusters = cs.get_connected_components(df.keys(),adj_list,dfSums)
    parentBarcodes,dfSums = cs.reduce_clusters(clusters,adj_list,dfSums)
    
    newDF = {}
    for parentBarcode,cluster in parentBarcodes.items():
        newSeries = pd.Series(0,index=genes,name=parentBarcode)
        for barcode in cluster:
            newSeries += df[barcode]
        if np.sum(newSeries) >= minCounts:
            newDF[parentBarcode] = newSeries

    print 'Collapsed_cells\t' + str(len(newDF))
    
    df = pd.DataFrame.from_dict(newDF)
    return df

# Write cell stats summary to file handle
def _writeCellStats(cell, cellbarcode, file_handle):
    file_handle.write('--------------------------------------------\n')
    file_handle.write('Cell barcode = ' + cellbarcode + '\n')
    file_handle.write('Number of cDNA transcripts = ' + str(cell.transcripts) + '\n')
    
    barcodeReads = sum(cell.gbcReadCounts.values())
    transcripts = sum(cell.gbcUMICounts.values())

    file_handle.write('Total gRNA Reads = ' + str(barcodeReads) + '\n')
    file_handle.write('Total gRNA transcripts = ' + str(transcripts) + '\n')
    
    file_handle.write('\nBarcodes by Transcript:\n')
    for barcode,count in sorted(cell.gbcUMICounts.items(),
                                key=lambda x: x[1], reverse=True):
        file_handle.write(barcode + '\t' + str(float(count)/transcripts) + '\n')
    
    file_handle.write('\nBarcodes by Read:\n') 
    for barcode,count in sorted(cell.gbcReadCounts.items(),
                                key=lambda x: x[1], reverse=True):
        file_handle.write(barcode + '\t' + str(float(count)/barcodeReads) + '\n')

# Print QC statistics on genotyped counts matrix
def qualityCheck(cell_barcodes_file, genotyped_df_file, gRNA_bam_file, cDNA_bam_file):
     
    total_cDNA_reads = int(str.strip(sp.check_output('samtools view -c ' + cDNA_bam_file,
                                                     shell = True)))
    mapped_cDNA_reads = int(str.strip(sp.check_output('samtools view -c -F 4 ' + cDNA_bam_file,
                                                      shell = True)))

    total_gRNA_reads = int(str.strip(sp.check_output('samtools view -c ' + gRNA_bam_file,
                                                     shell = True)))
    with open(cell_barcodes_file) as f:
        cell_barcodes = cp.load(f)
    
    df = pd.read_csv(genotyped_df_file, sep = '\t', header = 0, index_col = 0)
    
    transcripts = df.sum(0).to_dict()
    expressed_genes = df.astype(bool).sum(axis=0).to_dict()
    usable_cells = 0
    cell_grna_reads = {}
    cell_grna_transcripts = {}
    for barcode,cell in cell_barcodes.items():
        if cell.type == 'single':
            cell_grna_transcripts[barcode] = np.sum(cell.gbcUMICounts.values())
            cell_grna_reads[barcode] = np.sum(cell.gbcReadCounts.values())
            usable_cells += 1
    
    print 'Summary Stats'
    print 'Total cDNA reads:\t' + str(total_cDNA_reads)
    print 'cDNA mapping rate:\t' + str(float(mapped_cDNA_reads)/total_cDNA_reads)
    print 'Median cDNA UMIs per cell:\t' + str(np.median(transcripts.values()))
    print 'Median genes expressed per cell:\t' + str(np.median(expressed_genes.values()))
    print 'Total usable cells:\t' + str(usable_cells)
    print 'Total gRNA reads:\t' + str(total_gRNA_reads)
    print 'gRNA read usage rate:\t' + str(float(np.sum(cell_grna_reads.values()))/total_gRNA_reads)
    print 'Median gRNA UMIs per usable cell:\t' + str(np.median(cell_grna_transcripts.values()))
    print 'Median gRNA reads per usable cell:\t' + str(np.median(cell_grna_reads.values()))
    print 'gRNA clonal rate:\t' + str(1.0 - float(np.sum(cell_grna_transcripts.values()))/np.sum(cell_grna_reads.values()))
    

# Input: name of cell readcount file, number of cells to plot on x axis
# Ouput: knee plot for determining number of true cells in DropSeq
def findKnee(readCountFile1,outFileName,xlim,showPlot=True,
             readCountFile2=None,readCountFile3=None):
    vec1 = __parseReadCountFile__(readCountFile1, numParse = xlim, edit_dist=1)
    
    if readCountFile2:
        vec2 = __parseReadCountFile__(readCountFile2, numParse = xlim)
    if readCountFile2 and readCountFile3:
        vec3 = __parseReadCountFile__(readCountFile3, numParse = xlim)
    
    if not readCountFile2 and not readCountFile3:
        x = range(0,len(vec1))
        plt.plot(x,vec1)
    elif readCountFile2 and not readCountFile3:
        vecLen = min(len(vec1),len(vec2))
        x = range(0,vecLen)
        vec1 = vec1[0:vecLen]
        vec2 = vec2[0:vecLen]
        plt.plot(x,vec1,'b',x,vec2,'r')
    elif readCountFile2 and readCountFile3:
        vecLen = min(len(vec1),len(vec2),len(vec3))
        x = range(0,vecLen)
        vec1 = vec1[0:vecLen]
        vec2 = vec2[0:vecLen]
        vec3 = vec3[0:vecLen]
        plt.plot(x,vec1,'b',x,vec2,'r',x,vec3,'k')

    plt.xlabel('Cell Barcodes (sorted)')
    plt.ylabel('Cumulative Read Fraction')
    plt.xlim([0,len(vec1)])
    plt.savefig(outFileName)
    if showPlot:
        plt.show()

def __parseReadCountFile__(readCountFile,numParse,edit_dist=1):
    vec = pd.read_table(readCountFile,sep='\t',header=None,index_col=1)
    
    totalReads = float(np.sum(vec))
    
    cellCounts = {}
    for i in range(0,numParse):
        cellCounts[vec.index[i]] = int(vec.iloc[i])
    
    cellCounts = collapseBarcodesExact(cellCounts, edit_dist = edit_dist)
    print 'Number of cell barcodes = ' + str(len(cellCounts))

    vec = np.array(sorted(cellCounts.values(),reverse=True))
    vec = np.cumsum(vec)
    vec = np.divide(vec,totalReads)
    
    return vec


# Input: transcript matrix file in tsv format (must have 'gRNA' as gene)
# Output: histogram of gRNA transcripts per cell
def countGRNA(readCountDataFrameFile,countSummaryFile,nCore):
    transcriptCounts = {}
    with open(countSummaryFile) as f:
        f.readline()
        f.readline()
        f.readline()
        for line in f:
            fields = str.strip(line).split('\t')
            if len(fields) == 3:
                transcriptCounts[fields[0]] = int(fields[2])
    transcriptCounts = collapseBarcodesExact(transcriptCounts,1)
    transcriptCounts = sorted(transcriptCounts.items(),key = lambda x: x[1],reverse=True)
    transcriptCounts = transcriptCounts[0:nCore]
    coreBarcodes = [x[0] for x in transcriptCounts]

    df = pd.read_csv(readCountDataFrameFile,header=0,index_col=0,sep='\t')
    df = df[coreBarcodes]
    counts = df.loc['gRNA']
    
    
    print 'Cells: ' + str(len(counts))
    print 'Mean: ' + str(np.mean(counts))
    print 'Median: ' + str(np.median(counts))
    print 'Stdev: ' + str(np.std(counts))
    print 'Nonzero: ' + str(np.count_nonzero(counts))
    
    sns.distplot(counts,kde=False)
    plt.xlabel('gRNA Transcripts')
    plt.ylabel('Frequency')
    plt.show()

# Input: mouse count file, human count file, in tsv format;
#        number of core barcodes for each species
# Output: Species mixing plot
def speciesMixing(humanFile, mouseFile, nBarcodes = 50):
    mouseVec = pd.read_table(mouseFile, sep='\t', skiprows=2, header=0, index_col=0)
    humanVec = pd.read_table(humanFile, sep='\t', skiprows=2, header=0, index_col=0)
    
    mouseVec = mouseVec['NUM_TRANSCRIPTS']
    humanVec = humanVec['NUM_TRANSCRIPTS']
    
    combinedTable = pd.concat([humanVec, mouseVec], axis=1)
    combinedTable.columns = ['Human','Mouse']
    coreBarcodes = combinedTable.sum(1).sort(ascending=False,inplace=False).index[0:nBarcodes]
    combinedTable = combinedTable.loc[coreBarcodes]

    humanCells = {}
    mouseCells = {}
    doublets = 0
    colors = []
    for i,ind in enumerate(combinedTable.index):
        row = combinedTable.loc[ind]
        if row[0]/np.float(np.sum(row)) >= 0.95:
            humanCells[ind] = int(row[0])
            colors.append('r')
        elif row[1]/np.float(np.sum(row)) >= 0.95:
            mouseCells[ind] = int(row[1])
            colors.append('b')
        else:
            doublets += 1
            colors.append('k')

    numSingles = 'Percent Doublets: ' + str(2*float(doublets)/len(combinedTable.index)*100) + '%' 
    print numSingles
    print 'Number human = ' + str(len(humanCells))
    print 'Number mouse = ' + str(len(mouseCells))
    print 'Number of barcodes: ' + str(len(combinedTable.index))
    
    combinedTable = combinedTable.fillna(0)
    
    humanCells = sorted(humanCells.items(),key=lambda x: x[1],reverse=True)
    mouseCells = sorted(mouseCells.items(),key=lambda x: x[1],reverse=True)
    
    with open('humanTranscripts.txt','w') as f:
        for k,v in humanCells:
            f.write(k + '\t' + str(v) + '\n')
    
    with open('mouseTranscripts.txt','w') as f:
        for k,v in mouseCells:
            f.write(k + '\t' + str(v) + '\n')

    plt.scatter(combinedTable['Mouse'],combinedTable['Human'], c=colors)
    plt.xlabel('Mouse Transcripts')
    plt.ylabel('Human Transcripts')
    plt.title(numSingles) 
    plt.show()
 

if __name__ == '__main__': main()
