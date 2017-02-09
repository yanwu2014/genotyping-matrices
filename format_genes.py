import sys

genes_file = sys.argv[1]
output_file = genes_file.replace('.csv','.tsv')

f_out = open(output_file,'w')
with open(genes_file) as f:
    f.readline()
    for line in f:
        fields = str.strip(line).split(',')
        outline = fields[1].replace('_','\t')
        f_out.write(outline + '\n')

f_out.close()

