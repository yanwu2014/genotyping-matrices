import cPickle as cp
import numpy as np
import os
import multiprocessing as mp
from functools import partial
import time
import subprocess as sp
import pandas as pd
import sys
import shlex
from collections import defaultdict,Counter
import copy

def main():
    pass

# Annotated matrix of read counts per sample across a defined set of genes/transcripts
class CellMatrix:
        
    def __init__(self, countFileFolder, sampleFile = None, n_processors=1, 
                 level = 'gene', readLength=50, type='Salmon'):
        self.readLength = readLength

        # Get all count file names
        if sampleFile:
            countFileNames = []
            with open(sampleFile) as f:
                for line in f:
                    fields = str.strip(line).split('\t')
                    countFileNames.append(fields[0])
        else:
            countFileNames = os.listdir(countFileFolder)
            countFileNames = [x.replace('.counts.tsv','') for x in countFileNames]

        # Count reads for each cell
        t1 = time.time()
        pool = mp.Pool(n_processors)
        if str.lower(type) == 'salmon':
            _countFeatures_ = partial(countFeaturesSalmon,outerDirName=countFileFolder,level=level)
        else:
            _countFeatures_ = partial(countFeatures,countFileFolder=countFileFolder)
        countDFs = pool.map_async(_countFeatures_, countFileNames).get(9999999)
        t2 = time.time()
        print 'Read Counting Time = ' + str(t2-t1)

        # Get column and index names
        colNames = countFileNames
        rowNames = [str(x) for x in countDFs[0].index]

        # Create Pandas dataframe
        countMatrix = pd.DataFrame(index=rowNames,columns=colNames)
        for df in countDFs:
            countMatrix[df.columns.values[0]] = df[df.columns.values[0]]
        self.countMatrix = countMatrix

        # Calculate total reads
        totalReads = np.array(countMatrix.sum(0))
        self.totalReads = totalReads
        print 'Post-processing time = ' + str(time.time() - t2)
    


# Count features using count files generated from htseq-count script
def countFeatures(countFileName,countFileFolder):
    counts = {}
    countFilePath = countFileFolder + '/' + countFileName + '.counts.tsv'
    
    with open(countFilePath) as f:
        header = f.readline()
        for line in f:
            fields = line.split('\t')
            if str.strip(line) == '' or '__' in fields[0]:
                continue
            name = fields[0]
            if 'ENS' in name:
                gene = str.strip(name).split('.')[0]
            else:
                gene = str.strip(name)
            counts[gene] = int(fields[1])

    countsDF = pd.DataFrame.from_dict(counts,orient='index')
    countsDF.columns = [countFileName] 
    
    return countsDF 

# Count features using count files generated from salmon
def countFeaturesSalmon(countDirName, outerDirName, level = 'gene'):
    counts = Counter()
    countFilePath = outerDirName + '/' + countDirName + '/quant.sf'
    
    f = open(countFilePath)
    header = f.readline()
    for line in f:
        if line[0] == '#':
            continue
        fields = str.strip(line).split('\t')
        name = fields[0]
        tpm = float(fields[2])
        readCounts = float(fields[3])
        
        transcript = str.strip(name.split('|')[4])
        gene = str.strip(name.split('|')[5])
        
        if level == 'gene':
            counts[gene] += round(readCounts,0)
        elif level == 'transcript':
            counts[transcript] += round(readCounts,0)
        else:
            raise ValueError('Invalid count level')

    f.close()
    
    countsDF = pd.DataFrame.from_dict(counts,orient='index')
    countsDF.columns = [countDirName] 
    
    return countsDF
        
if __name__=='__main__':main()
