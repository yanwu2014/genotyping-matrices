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

df = pd.read_csv(dataFrameFile, sep = '\t', header = 0, index_col = 0)
coreCells = df.sum(0).to_dict()

cellBarcodes = sc.getCellBarcodes(bamFile, coreCells, edit_dist = cbc_edit_dist)
for cbc,cell in cellBarcodes.items():
    cellBarcodes[cbc].collapseMolTags(cbc_edit_dist)

barcodeFile = bamFile.replace('.sorted.tagged.bam', '_cell_barcodes.pickle')
with open(barcodeFile, 'w') as f:
    cp.dump(cellBarcodes,f)


