### Wrapper script for extracting gRNA reads ###

import sys
import cPickle as cp
import spacerCalling as sc
import pandas as pd


dataFrameFile = sys.argv[1]
bamFile = sys.argv[2]
BARCODE_LENGTH = int(sys.argv[3])
BC_START_HANDLE = sys.argv[4]
BC_END_HANDLE = sys.argv[5]

cbc_edit_dist = 1 # Edit distance between cells
gbc_edit_dist = 1

#BARCODE_LENGTH = 20

## gRNA construct handles
#BC_START_HANDLE = 'CCGAGTCGGTGC'
#BC_END_HANDLE = 'TATGA'

## ORF overexpression handles
#BC_START_HANDLE = 'GGCTGTTACGCG'
#BC_END_HANDLE = 'CTACTGAC'

# Load scRNA-seq counts matrix
df = pd.read_csv(dataFrameFile, sep = '\t', header = 0, index_col = 0)

# Strip batch ids
cell_barcodes = [x.split("-")[0] for x in df.columns]
df.columns = cell_barcodes

# Get cell names and total UMIs
coreCells = df.sum(0).to_dict()

# Parse cell barcodes
cellBarcodes = sc.getCellBarcodes(bamFile, coreCells, BARCODE_LENGTH, BC_START_HANDLE,
                                  BC_END_HANDLE, edit_dist = cbc_edit_dist)

# Collapse barcodes within edit distance
for cbc,cell in cellBarcodes.items():
    cellBarcodes[cbc].collapseMolTags(cbc_edit_dist)

# Write barcodes dictionary to pickled object
barcodeFile = bamFile.replace('.tagged.bam', '_cell_barcodes.pickle')
with open(barcodeFile, 'w') as f:
    cp.dump(cellBarcodes,f)

