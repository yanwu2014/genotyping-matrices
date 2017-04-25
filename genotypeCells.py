# This is a wrapper script for genotyping cells
import sys
import pandas as pd
from dropSeqModule import *
import spacerCalling as sc

dataFrameFile = sys.argv[1]
barcodeFile = sys.argv[2] 
barcodeToGeneFile = sys.argv[3]

readThreshFrac = float(sys.argv[4])
umiThreshFrac = float(sys.argv[5])

minUMIFrac = float(sys.argv[6])

with open(barcodeFile) as f:
    cellBarcodes = cp.load(f)

cellBarcodes = callBarcodes(dataFrameFile, cellBarcodes, barcodeToGeneFile, readThreshFrac,
                            umiThreshFrac, minUMIFrac = minUMIFrac)

