### Wrapper script for mapping barcode to spacer from plasmid sequencing ###

import sys
from collections import defaultdict,Counter
import spacerCalling as sc
from time import time
import cPickle as cp

def main():
    plasmidFastqFile = sys.argv[1]
    plasmidFastqFileMate = sys.argv[2]
    spacerFile = sys.argv[3]
    edit_dist = int(sys.argv[4])
    read_threshold = int(sys.argv[5])

    t1 = time()
    barcodeToSpacer = sc.parsePlasmidFile(plasmidFastqFile,spacerFile,edit_dist,
                                          plasmidFastqFileMate,read_threshold)
    print 'Time Elapsed = ' + str(time() - t1)
     
    name = sys.argv[1].split('/')[-1].split('_')[0]
    with open(name + '_barcode_to_gene.pickle','w') as f:
        cp.dump(barcodeToSpacer,f)

if __name__ == '__main__': main()
