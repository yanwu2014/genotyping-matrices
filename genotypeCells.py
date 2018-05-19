# This is a wrapper script for genotyping cells
import sys
import pandas as pd
from spacerCalling import *

barcodeFile = sys.argv[1]
outputFile = sys.argv[2]
barcodeToGeneFile = sys.argv[3]

readThreshFrac = float(sys.argv[4])
umiThreshFrac = float(sys.argv[5])

minFrac = float(sys.argv[6])

with open(barcodeFile) as f:
    cellBarcodes = cp.load(f)

barcodeToGene = {}
with open(barcodeToGeneFile) as f:
    for line in f:
        fields = str.strip(line).split("\t")
        barcodeToGene[fields[0]] = fields[1]

if minFrac >= 1.0:
    minReads = int(minFrac)
    genotype_dict = callBarcodes(cellBarcodes, barcodeToGene, readThreshFrac,
                                 umiThreshFrac, minUMIReads = minReads)
else:
    genotype_dict = callBarcodes(cellBarcodes, barcodeToGene, readThreshFrac,
                                 umiThreshFrac, minUMIFrac = minFrac)

with open(outputFile, 'w') as f:
    for genotype,cells in genotype_dict.items():
        outLine = genotype + ',\"' + ",".join(cells) + '\"' + '\n'
        f.write(outLine)

noPlasmids = 0
gbc_dist = Counter()
for cell in cellBarcodes.values():
    if cell.type == 'usable':
        gbc_dist[len(cell.genotype)] += 1
    elif cell.type == 'noGRNA':
        gbc_dist[0] += 1
    else:
        noPlasmids += 1

print 'Total_Cells\t' + str(len(cellBarcodes))
print 'No_Plasmid\t' + str(noPlasmids)
for k,v in gbc_dist.items():
    print str(k) + "\t" + str(v)
