### Script to rename pickled python object ###

import cPickle as cp
import sys
import subprocess as sp

with open(sys.argv[1]) as f:
    object = cp.load(f)

with open(sys.argv[2],'w') as f:
    cp.dump(object,f)

sp.call('rm ' + sys.argv[1],shell=True)
