import sys
import subprocess as sp
import multiprocessing as mp
import os
from functools import partial

def main():
    targetFolder = sys.argv[1]
    outFolder = sys.argv[2]
    sampleFile = sys.argv[3]
    annotationFile = sys.argv[4]
    n_processors = int(sys.argv[5])
    
    if str.upper(sampleFile) == 'NONE':
        fileNames = str.strip(sp.check_output('ls ' + targetFolder + '/*.bam',shell=True)).split('\n')
        fileNames = [x.split('/')[-1].replace('.Aligned.sortedByCoord.out.bam','') for x in fileNames]
    
    else:
        fileNames = []
        with open(sampleFile) as f:
            for line in f:
                fields = str.strip(line).split('\t')
                fileNames.append(fields[0])
    
    partialRunHTseq = partial(runHTseq,targetFolder=targetFolder, outFolder=outFolder, 
                              annotationFile=annotationFile)

    pool = mp.Pool(n_processors)
    pool.map_async(partialRunHTseq,fileNames).get(9999999)


def runHTseq(fileName, targetFolder, outFolder, annotationFile):
    filePath = targetFolder + '/' + fileName + '.Aligned.sortedByCoord.out.bam'
    outFilePath = outFolder + '/' + fileName + '.counts.tsv'
    cmd = 'htseq-count -a 20 -f bam -s no -t exon ' + filePath + ' ' + \
          annotationFile + ' > ' + outFilePath
    sp.call(cmd, shell=True)

if __name__=='__main__': main() 
