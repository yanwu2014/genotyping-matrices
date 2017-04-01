import sys
import subprocess as sp
import os
from functools import partial

def main():
    targetFolder = sys.argv[1]
    outFolder = sys.argv[2]
    sampleFile = sys.argv[3]
    indexFile = sys.argv[4]
    n_processors = int(sys.argv[5])
     
    if str.upper(sampleFile) == 'NONE':
        fileNames = str.strip(sp.check_output('ls ' + targetFolder + '/*.fastq.gz',shell=True)).split('\n')
        fileNames = [x.split('/')[-1].replace('.fastq.gz','') for x in fileNames]
    
    else:
        fileNames = []
        with open(sampleFile) as f:
            for line in f:
                fields = str.strip(line).split('\t')
                fileNames.append(fields[0])
    
    for file in fileNames:
        runSalmon(file,targetFolder,outFolder,indexFile) 

def runSalmon(fileName, targetFolder, outFolder, indexFile):
    filePath = targetFolder + '/' + fileName + '.fastq.gz'
    outDirPath = outFolder + '/' + fileName
    cmd = "salmon quant -i " + indexFile + " -l U -r " + "<(zcat " + filePath + ") -o " + outDirPath
    sp.call(['/bin/bash','-c',cmd])

if __name__=='__main__': main() 
