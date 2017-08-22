### Wrapper script for extracting gRNA reads ###

import sys
import cPickle as cp
import spacerCalling as sc
import pandas as pd
#import dropSeqModule as ds


dataFrameFile = sys.argv[1]
bamFile = sys.argv[2]
cbc_edit_dist = int(sys.argv[3]) # Edit distance between cells
gbc_edit_dist = 1

BARCODE_LENGTH = 20
#BC_START_HANDLE = 'CCGAGTCGGTGC'
#BC_END_HANDLE = 'TATGA' 
BC_START_HANDLE = 'GGCTGTTACGCG'
BC_END_HANDLE = 'CTACTGAC'


df = pd.read_csv(dataFrameFile, sep = '\t', header = 0, index_col = 0)

# Strip batch ids
cell_barcodes = [x.split("-")[0] for x in df.columns]
df.columns = cell_barcodes

# Get cell names and total UMIs
coreCells = df.sum(0).to_dict()

cellBarcodes = sc.getCellBarcodes(bamFile, coreCells, BARCODE_LENGTH, BC_START_HANDLE,
                                  BC_END_HANDLE, edit_dist = cbc_edit_dist)
for cbc,cell in cellBarcodes.items():
    cellBarcodes[cbc].collapseMolTags(cbc_edit_dist)

barcodeFile = bamFile.replace('.tagged.bam', '_cell_barcodes.pickle')
with open(barcodeFile, 'w') as f:
    cp.dump(cellBarcodes,f)


