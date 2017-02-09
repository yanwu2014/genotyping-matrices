### Wrapper script for running quality check on genotyped dataframe ###

import sys
import dropSeqModule as ds

# Arguments
# 1: cell_barcodes pickle file
# 2: genotyped dataframe
# 3: gRNA bam file
# 4: cDNA bam file

ds.qualityCheck(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
