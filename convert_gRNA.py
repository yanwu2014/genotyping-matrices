import sys
import pysam as ps
import subprocess as sp
import pandas as pd
from dedupBarcodes import exactExists,belowQual
from collections import Counter

## Script to convert 10X gRNA fastq files to a tagged bam file
## Each gRNA read in the aligned bam file has a cell/molecule tag
def main():
    read1_fastq_file = sys.argv[1]
    read2_bam_file = sys.argv[2]
    tag_reads(read1_fastq_file, read2_bam_file)
    
# Tags bam file with cell/molecule barcodes from read 1 fastq file
# Read 1: cell barcode + UMI
# Read 2: cDNA
def tag_reads(read1_fastq_file, read2_bam_file):
    cellTags = get_tags(read1_fastq_file, 0, 16)
    molTags = get_tags(read1_fastq_file, 16, 26)
    
    f_in = ps.AlignmentFile(read2_bam_file, 'rb')
    f_out = ps.AlignmentFile(read2_bam_file.replace('.bam','.tagged.bam'), 'wb',
                             template = f_in)
    for read in f_in.fetch():
        if read.query_name in cellTags and read.query_name in molTags:
            read.set_tag('XM',molTags[read.query_name])
            read.set_tag('XC',cellTags[read.query_name])
            f_out.write(read)

    f_in.close()
    f_out.close()

# Parse read file and extract sequence from start to end
def get_tags(read_file, start, end):
    celltags = {}
    fh = ps.FastxFile(read_file)
    for read in fh:
        tag = read.sequence[start:end]
        quals = list(read.quality)[start:end]
        quals = map(lambda x: ord(x) - 33, quals)
        if belowQual(quals,10) < 2:
            celltags[read.name] = tag
    fh.close()
    return celltags

# Combine gzipped fastq files
def cat_files(sample_ids, wildcard, out_prefix):
    in_files = []
    for sample_dir in sample_ids:
        cmd = 'ls ' + sample_dir + wildcard
        in_paths = sp.check_output(cmd, shell = True)
        in_paths = str.strip(in_paths).split('\n')

        in_files.extend(in_paths)

    cat_str = ''
    for fi in in_files:
        cat_str = fi + ' ' + cat_str
    cat_str = 'cat ' + cat_str + '> ' + out_prefix + '.fastq.gz'
    
    sp.call(cat_str, shell = True)

        
if __name__ == '__main__': main()
