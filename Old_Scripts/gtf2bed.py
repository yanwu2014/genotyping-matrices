import sys

# Input arguments
input_gtf = sys.argv[1]
output_bed = input_gtf.replace('.gtf','.bed')

with open(input_gtf) as f:
    for line in f:
        if line[0] == '#':
            continue
        fields = str.strip(line).split('\t')
        entryType = fields[2]
        
        if entryType == 'gene':
            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            
            metadata = fields[8].split(';')
            name = metadata[3].split(' ')[-1][1:-1]
            print chrom + '\t' + start + '\t' + end + '\t' + name
