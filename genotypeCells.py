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
minUMIreads = int(sys.argv[6])

with open(barcodeFile) as f:
    cellBarcodes = cp.load(f)

cellBarcodes = callBarcodes(dataFrameFile, cellBarcodes, barcodeToGeneFile,
                            edit_dist = edit_dist, singleFrac = singleFrac,
                            minUMIreads = minUMIreads)

