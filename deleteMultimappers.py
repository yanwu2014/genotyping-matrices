import sys
from collections import defaultdict

# Input: human mapped sam file, mouse mapped sam file
# Output: human & mouse sam files with reads mapping to both genomes removed
def trimMultimappers(hSamFile,mSamFile):
    hReads = getReadsSam(hSamFile)
    mReads = getReadsSam(mSamFile)
    
    multimappers = hReads & mReads 
    
    trimSam(hSamFile,multimappers)
    trimSam(mSamFile,multimappers)


def trimSam(samFile,multimappers):
    samFileOut = samFile.replace('.sam','.species_trimmed.sam')
    f_out = open(samFileOut,'w')
    with open(samFile) as f:
        for line in f:
            line = str.strip(line)
            if line[0] == '@':
                f_out.write(line + '\n')
            elif not str.strip(line.split('\t')[0]) in multimappers:
                f_out.write(line + '\n')
    f_out.close()


def getReadsSam(samFile):
    reads = set() 
    
    with open(samFile) as f:
        for line in f:
            if line[0] != '@':
                line = str.strip(line)
                fields = line.split('\t')
                id = str.strip(fields[0])
                reads.add(id)
    return reads

if __name__ == '__main__': main()
