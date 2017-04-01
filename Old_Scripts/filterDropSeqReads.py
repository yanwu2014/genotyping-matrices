from Bio import SeqIO
import sys
import gzip as gz

def main():
    fastqFile = sys.argv[1]
    fastqMate = sys.argv[2]

    handle = gz.open(fastqFile, "rU")
    mateHandle = gz.open(fastqMate,'rU')
    
    readsToKeep = set()
    discard = 0
    reads = 0
    for record in SeqIO.parse(handle, "fastq"):
        reads += 1
        if checkPolyT(str(record.seq),20,30) < 0.7:
            discard += 1
        else:
            readsToKeep.add(record.id)
    
    print 'Percent Reads Discarded = ' + str(float(discard)/reads)
    
    handle.seek(0,0)
    mateHandle.seek(0,0)

    trimmedRead1Iterator = (record for record in SeqIO.parse(handle,'fastq') \
                            if record.id in readsToKeep)
    trimmedRead2Iterator = (record for record in SeqIO.parse(mateHandle,'fastq') \
                            if record.id in readsToKeep)
    
    newRead1File = fastqFile.replace('.fastq.gz','.filtered.fastq.gz')
    newRead2File = fastqMate.replace('.fastq.gz','.filtered.fastq.gz')

    read1Output = gz.open(newRead1File,'w')
    read2Output = gz.open(newRead2File,'w')

    SeqIO.write(trimmedRead1Iterator,read1Output,'fastq')
    SeqIO.write(trimmedRead2Iterator,read2Output,'fastq')

    handle.close()
    mateHandle.close()
    read1Output.close()
    read2Output.close()


def checkPolyT(string,start,end):
    substring = string[start:end]
    return float(sum(map(lambda x: x == 'T',substring)))/len(substring)  

if __name__ == '__main__': main()
