# This is a wrapper script for genotyping cells
import sys
import pandas as pd
from spacerCalling import *

dataFrameFile = sys.argv[1]
barcodeFile = sys.argv[2] 
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
    cellBarcodes = callBarcodes(dataFrameFile, cellBarcodes, barcodeToGene, readThreshFrac,
                                umiThreshFrac, minUMIReads = minReads)
else:
    cellBarcodes = callBarcodes(dataFrameFile, cellBarcodes, barcodeToGene, readThreshFrac,
                                umiThreshFrac, minUMIFrac = minFrac)


