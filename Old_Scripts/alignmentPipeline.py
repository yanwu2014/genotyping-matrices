import sys
import os
import subprocess as sp
import shlex
import signal

def main():
    genomeDir = sys.argv[1] 
    folder = sys.argv[2]
    outFolder = sys.argv[3]
    nThreads = int(sys.argv[4])
    
    # Whether or not the reads are paired end
    pairedEnd = sys.argv[5]
    if pairedEnd == 'True':
        pairedEnd = True
    else:
        pairedEnd = False
    
    if pairedEnd:
        read1Files = str.strip(sp.check_output('ls ' + folder + '/*R1*',shell=True)).split('\n')
        read2Files = [x.replace('R1','R2') for x in read1Files]
        readMates = zip(read1Files,read2Files)
        
        for mates in readMates:
            runSTARpaired(mates, folder, genomeDir, outFolder, nThreads)
    else:
        files = os.listdir(folder)
        for file in files:
            runSTARsingle(file, folder, genomeDir, outFolder, nThreads)


def runSTARpaired(mates, folder, genomeDir, outFolder, nThreads):
    outFileName = outFolder + '/' + mates[0].split('/')[-1].replace('_R1','').replace('fastq.gz','')
    starCmd = 'STAR  --runThreadN ' + str(nThreads) + ' --genomeDir ' + \
          genomeDir + ' --readFilesIn ' + mates[0] + ' ' + mates[1] + ' --outSAMtype ' + \
          ' BAM SortedByCoordinate --outFileNamePrefix ' + outFileName + \
          ' --readFilesCommand zcat --genomeLoad LoadAndKeep ' + \
          ' --limitBAMsortRAM 30000000000 --limitIObufferSize 350000000' 
    sp.call(starCmd,shell=True)

def runSTARsingle(file, folder, genomeDir, outFolder, nThreads):
    filepath = folder + '/' + file 
    outFileName = outFolder + '/' + file.replace('fastq.gz','')
    starCmd = 'STAR --outSAMunmapped Within --runThreadN ' + str(nThreads) + ' --genomeDir ' + \
          genomeDir + ' --readFilesIn ' + filepath + ' --outSAMtype' + \
          ' BAM SortedByCoordinate --outFileNamePrefix ' + \
          outFileName + ' --outSAMstrandField intronMotif' + \
          ' --readFilesCommand zcat --genomeLoad LoadAndKeep ' + \
          ' --limitBAMsortRAM 30000000000 --limitIObufferSize 350000000'
    sp.call(starCmd,shell=True)

if __name__ == '__main__': main()
