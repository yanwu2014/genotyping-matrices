import sys
import pysam as ps
import numpy as np
import pandas as pd
from collections import defaultdict,Counter
from dedupBarcodes import exactExists,collapseBarcodesExact
from itertools import islice
import cPickle as cp

def main():
    bamFileName = sys.argv[1]
    fixedCutoff = float(sys.argv[2])
    minReads = int(sys.argv[3])

    totalReads = 0
    baseDist = defaultdict(Counter)
    bamFile = ps.AlignmentFile(bamFileName, 'rb')
    for read in bamFile.fetch():
        totalReads += 1

        molTag = read.get_tag('XM')
        cellTag = read.get_tag('XC')
        
        baseDist[cellTag][molTag[-1]] += 1

    bamFile.close()

    baseDistCounts = {k: sum(v.values()) for k,v in baseDist.items() if sum(v.values()) > minReads}
    baseDist = {k:v for k,v in baseDist.items() if k in baseDistCounts}
    
    #baseDist,baseDistCounts = collapseBarcodesExact(baseDist, edit_dist = 2, barcodeCounts = barcodeDistCounts)
    
    for k,counts in baseDist.items():
        baseDist[k] = {base:float(count)/sum(counts.values()) for base,count in counts.items()}

    needsCorrection = set()
    for cellTag,freqCounter in baseDist.items():
        if 'T' in freqCounter and freqCounter['T'] > fixedCutoff:
            needsCorrection.add(cellTag[0:11])

    readsCorrected = 0

    bamFile = ps.AlignmentFile(bamFileName, 'rb')
    outFile = ps.AlignmentFile(bamFileName.replace('.bam','.cleaned.bam'), 'wb', template = bamFile)
    for read in bamFile.fetch():
        molTag = read.get_tag('XM')
        cellTag = read.get_tag('XC')
        
        if cellTag[0:11] in needsCorrection:
            cellTag = list(cellTag)
            molTag = list(molTag)
            
            molTag[1:] = molTag[:-1]
            molTag[0] = cellTag[-1]
            molTag = ''.join(molTag)

            cellTag[-1] = 'N'
            cellTag = ''.join(cellTag)

            read.set_tag('XM',molTag)
            read.set_tag('XC',cellTag)
            
            outFile.write(read)
            readsCorrected += 1

        else:
            outFile.write(read)

    bamFile.close()
    outFile.close()
    
    print 'Total Cells = ' + str(len(baseDist))
    print 'Fraction Cells Needing Correction = ' + str(float(len(needsCorrection))/len(baseDist))
    print 'Fraction Reads Needing Correction = ' + str(float(readsCorrected)/totalReads)
    

# Input: data frame, desired edit distance to be collapsed
# Output: collapsed data frame where cell barcodes within edit_dist have their count arrays added 
def _collapseCellBarcodes(baseDist, baseDistCounts, edit_dist, minCounts = 200, discardError = True):
    cs = ClusterAndReducer()

    adj_list = cs.get_adj_list(df.keys(), dfSums, 2, isCellBarcode = False,
                               countRatio = 1.0)
    clusters = cs.get_connected_components(df.keys(),adj_list,dfSums)
    parentBarcodes,dfSums = cs.reduce_clusters(clusters,adj_list,dfSums)
    
    newDF = {}
    for parentBarcode,cluster in parentBarcodes.items():
        if parentBarcode[11] == 'N' and discardError:
            if np.sum(df[parentBarcode]) >= minCounts:
                newDF[parentBarcode] = df[parentBarcode]
        else:
            newSeries = pd.Series(0,index=genes,name=parentBarcode)
            for barcode in cluster:
                newSeries += df[barcode]
            if np.sum(newSeries) >= minCounts:
                newDF[parentBarcode] = newSeries

    df = pd.DataFrame.from_dict(newDF)
    return df


if __name__ == '__main__': main()

