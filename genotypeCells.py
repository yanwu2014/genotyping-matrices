# This is a wrapper script for genotyping cells
import sys
import pandas as pd
from dropSeqModule import *
import spacerCalling as sc

dataFrameFile = sys.argv[1]
barcodeFile = sys.argv[2] 
barcodeToGeneFile = sys.argv[3]
edit_dist = int(sys.argv[4])
singleFrac = float(sys.argv[5])
dualFrac = float(sys.argv[6])
tripFrac = float(sys.argv[7])
minUMIreads = int(sys.argv[8])

minGBCreads = 20
minGBCumis = 2

with open(barcodeFile) as f:
    cellBarcodes = cp.load(f)

cellBarcodes = callBarcodes(dataFrameFile, cellBarcodes, barcodeToGeneFile,
                            edit_dist = edit_dist, singleFrac = singleFrac,
                            dualFrac = dualFrac, tripFrac = tripFrac, 
                            minUMIreads = minUMIreads, minReads = minGBCreads,
                            minUMIs = minGBCumis)

