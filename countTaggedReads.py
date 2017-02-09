### Script to count how many reads map to Introns, Exons, and Intergenic Regions

import pysam as ps
import sys
from itertools import islice

introns = 0
exons = 0
intergenic = 0
totalReads = 0
samFile = ps.AlignmentFile(sys.argv[1])
for read in samFile.fetch(until_eof=True):
    totalReads += 1
    if read.has_tag('XF'):
        mappedLoc = read.get_tag('XF')
        if mappedLoc == 'INTRONIC':
            introns += 1
        elif mappedLoc == 'INTERGENIC':
            intergenic += 1
        elif mappedLoc == 'CODING':
            exons += 1


print 'Total Reads = ' + str(totalReads)
print 'Intronic Reads = ' + str(introns)
print 'Exonic Reads = ' + str(exons)
print 'Intergenic Reads = ' + str(intergenic)
