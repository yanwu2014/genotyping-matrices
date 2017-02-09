# Wrapper script for making DropSeq knee plots
import matplotlib
matplotlib.use('Agg')
import dropSeqModule as ds
import sys

readCountFile = sys.argv[1]
outFileName = sys.argv[2]
numCells = int(sys.argv[3])

ds.findKnee(readCountFile,outFileName,4*numCells,showPlot=False)
