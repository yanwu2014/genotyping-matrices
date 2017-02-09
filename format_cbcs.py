import sys

barcodes_file = sys.argv[1]
output_file = barcodes_file.replace('.csv','.tsv')

f_out = open(output_file,'w')
with open(barcodes_file) as f:
    f.readline()
    for line in f:
        fields = str.strip(line).split(',')
        f_out.write(fields[1] + '\n')
f_out.close()

