### Primary module for genotyping cells ###

import sys
import cPickle as cp
from editdistance import eval
from collections import defaultdict,Counter
from random import sample
from Bio import SeqIO
from hammingSearch import approxHammingSearch 
from random import random
from time import time
import gzip
import pysam as ps
from itertools import islice
from operator import itemgetter
from dedupBarcodes import *
import distance
import numpy as np
import pandas as pd

BARCODE_LENGTH = 12
#BC_START_HANDLE = 'CCGAGTCGGTGC'
#BC_END_HANDLE = 'TATGA' 
BC_START_HANDLE = 'GGCTGTTACGCG'
BC_END_HANDLE = ''

cs = ClusterAndReducer()

# Cell class
class Cell:
    
    # Construct Cell object with cell barcode
    def __init__(self, cellTag):
        self.transcripts = 0 # Number of cDNA UMIs
        self.gbc_reads = 0 # Number of gRNA reads
        
        self.type = None # Type: no-Plasmid, no-gRNA, usable
        self.genotype = None # Genotype: e.g. HDAC2-1,HSP90AA1-1
        
        self.cellTag = cellTag # Cell Barcode
        self.molTags = defaultdict(Counter) # Barcode - UMI - ReadCount dict
        self.gbcReadCounts = Counter() # Barcode - ReadCount dict
        self.gbcUMICounts = Counter() # Barcode - UMIcount dict

    # Input: molecule tag, corresponding guide barcode
    def addMolTag(self, molTag, guideBarcode):
        self.molTags[guideBarcode][molTag] += 1
        self.gbc_reads += 1
    
    # Collapse molecule and UMI tags within edit distance
    def collapseMolTags(self, edit_dist):
        guideTags = self.molTags
        
        for gbc, molTagCounts in guideTags.items():
            guideTags[gbc] = Counter(collapseBarcodesExact(molTagCounts, edit_dist=edit_dist,
                                                           hamming = True))
        
        guideTagCounts = _sumDict(guideTags) 
        guideTags,guideTagCounts = collapseBarcodesExact(guideTags, edit_dist, 
                                                         barcodeCounts = guideTagCounts)
           
        self.molTags = guideTags


    # Input: self.molTags, edit distance, minReadFrac
    # Output: count guide barcodes and store in barcodeCounts
    def countBarcodes(self, minReads = None, minFrac = None):
        guideTags = self.molTags
        if minReads:
            guideTags = _thresholdBarcodes(guideTags, minReads)
        elif minFrac:
            minReads = round(minFrac*self.gbc_reads)
            guideTags = _thresholdBarcodes(guideTags, minReads)

        for gbc, molTagCounts in guideTags.items(): 
            self.gbcReadCounts[gbc] = sum(molTagCounts.values())
            self.gbcUMICounts[gbc] = len(molTagCounts)
        
        self.molTags = guideTags
    
    # Takes the read/umi counts calculated from countBarcodes, the barcode to gene
    # dictionary, and attempts to genotype each cell with a minimum read and UMI
    # fraction threshold
    def callBarcode(self, barcodeToGene, readThresh, umiThresh, minReads = 5, minUMIs = 2):
        # Map each barcode to a guide RNA and count
        guideReadCounts = Counter()
        guideUMICounts = Counter()
        for gbc in self.gbcReadCounts:
            libGBC = exactExists(barcodeToGene, gbc, edit_dist = 1)
            if libGBC and type(libGBC) == str:
                guide = barcodeToGene[libGBC]
                guideReadCounts[guide] += self.gbcReadCounts[gbc]
                guideUMICounts[guide] += self.gbcUMICounts[gbc]
        
        if len(guideReadCounts) == 0:
            self.type = "noPlasmid"
            self.genotype = ""
        
        genotypes = []
        totalReads = np.sum(guideReadCounts.values())
        totalUMIs = np.sum(guideUMICounts.values())
        for guide in guideReadCounts:
            reads = guideReadCounts[guide]
            umis = guideUMICounts[guide]
            if float(reads)/totalReads > readThresh and float(umis)/totalUMIs > umiThresh and \
                    reads >= minReads and umis >= minUMIs:
                genotypes.append(guide)
         
        if len(genotypes) == 0:
            self.type = 'noGRNA'
            self.genotype = ''

        else:
            self.type = 'usable'
            self.genotype = genotypes

            
# Input: mapped & corrected BAM file from DropSeq, barcodeToSpacer mapping dict, 
#        cell readcount file from DropSeq, number of desired core barcodes
# Output: dict with cell barcodes as keys, Cell objects as vals 
def getCellBarcodes(bamFile, coreCells, edit_dist = 1):
    
    bamFile = ps.AlignmentFile(bamFile,'rb') # load bam file
    
    cellBarcodes = {} 
    nreads = 0
    for read in bamFile.fetch('HUMAN_MOUSE_chrGRNA'):
        molTag = _getTag(read,'XM')
        cellTag = _getTag(read,'XC')
        
        guideBarcode = _getBarcode(read)
        if guideBarcode:
            nreads += 1
            if cellTag in cellBarcodes:
                cellBarcodes[cellTag].addMolTag(molTag, guideBarcode)
            else:
                cell = Cell(cellTag)
                cell.addMolTag(molTag, guideBarcode)
                cellBarcodes[cellTag] = cell 
            
            if (nreads % 1000000) == 0: print 'Reads parsed = ' + str(nreads)
    
    bamFile.close()
    print 'gRNA Reads_Used\t' + str(nreads)
    
    # Match cell barcodes to known cell barcodes in counts matrix
    cellBarcodes = _collapseCells(cellBarcodes, coreCells, edit_dist = edit_dist)
        
    return cellBarcodes

# Maps the cell barcodes from the gRNA calling onto the cells from the counts matrix
def _collapseCells(cellBarcodes, coreCells, edit_dist):
    cellTagCollisions = 0
    doesExist = 0
    cellsToCollapse = defaultdict(set)
    print 'Original_Cell_Tags\t' + str(len(cellBarcodes))
    for cellTag in cellBarcodes:
        cellTagExists = exactExists(coreCells, cellTag, edit_dist = edit_dist)

        if cellTagExists:
            if type(cellTagExists) == str:
                cellsToCollapse[cellTagExists].add(cellTag)
                doesExist += 1
            else:
                cellTagCollisions += 1
                #bestCellTag = argmax([(tag,coreCells[tag]) for tag in cellTagExists],
                #                     scorePos = 1)
                #cellsToCollapse[bestCellTag].add(cellTag)

    print 'Cell_Tag_Collisions\t' + str(cellTagCollisions)
    print 'Cells_matched_to_core_cells\t' + str(doesExist)
    
    collapsedCellBarcodes = {}
    for parentBC,children in cellsToCollapse.items():
        newCell = _mergeCells(cellBarcodes, parentBC, children)
        newCell.transcripts = coreCells[parentBC] 
        collapsedCellBarcodes[parentBC] = newCell

    return collapsedCellBarcodes

# Merge the cells mapping to the parent barcode (from the counts matrix)
def _mergeCells(cellBarcodes, parentBC, children):
    newCell = Cell(parentBC)
    for childBC in children:
        childCell = cellBarcodes[childBC]
        newCell.molTags = _mergeCounterDict(newCell.molTags,childCell.molTags)
        newCell.gbc_reads = newCell.gbc_reads + childCell.gbc_reads
    return newCell

# Merge two dictionaries of counters (or counters)
def _mergeCounterDict(counterDictA, counterDictB):
    for key,counter in counterDictB.items():
        if key in counterDictA:
            counterDictA[key] += counter
        else:
            counterDictA[key] = counter
    return counterDictA


# Input: pysam AlignedSegment object, SAM tagName
# Output: SAM tag if it exists, False otherwise 
def _getTag(read,tagName):
    for tag,val in read.tags:
        if tag == tagName: 
            return val
    return False
                    
# Input: pysam AlignedSegment object
# Output: gRNA barcode if it exists, False otherwise
def _getBarcode(read):
    # Make sure read maps to end of reference (right before barcode)
    if len(read.query_sequence) < 65:
        return False
    
    # ensure that the sequence right before the barcode is what we expect
    bcStart = list(approxHammingSearch(BC_START_HANDLE, read.query_sequence))
    if len(bcStart) < 1: 
        return False
    left_pointer = argmin(bcStart) + len(BC_START_HANDLE)
    
    if len(BC_END_HANDLE) > 1:
        bcEnd = list(approxHammingSearch(BC_END_HANDLE, read.query_sequence))
        if len(bcEnd) < 1:
            return False
        right_pointer = argmin(bcEnd)
    else:
        right_pointer = left_pointer + BARCODE_LENGTH

    barcode = read.query_sequence[left_pointer:right_pointer]

    # Ensure the read covers the entire barcode
    if (len(barcode) < BARCODE_LENGTH - 1) or (len(barcode) > BARCODE_LENGTH):
        return False

    quals = read.query_qualities[left_pointer:right_pointer]
    if belowQual(quals,20) > 1:
        return False

    return barcode

# Threshold counter dict by minimum ratio of total counts
def _thresholdCounterFrac(counter, ratio = 0.01):
    totalCount = sum(counter.values())
    for k,count in counter.items():
        if float(count)/totalCount < ratio:
            del counter[k]
    return counter

# Threshold counter dict by minimum counts
def _thresholdCounterRead(counter, minReads = 1):
    for k,count in counter.items():
        if count < minReads:
            del counter[k]
    return counter

### Modules for mapping barcodes to spacers from plasmid amplicon sequencing ###

# Input: Name of CRISPR plasmid library Fastq file linking spacers to gRNA barcodes
# Output: Dict with barcode as key and spacer as value. 
def parsePlasmidFile(plasmidFastqFile,spacer_file,edit_dist=1,plasmidFastqFileMate=None,
                     threshold=2):
    spHandles = ('ACGAAACACCG','GTTTTAGAGC')
    bcHandles = ('CGAGTCGGTGC','TATGAGGA')

    if plasmidFastqFileMate:
        barcodeSpacerCounts = _parsePairedEnd(plasmidFastqFile,plasmidFastqFileMate,
                                              spHandles,bcHandles,edit_dist)
    else:
        barcodeSpacerCounts = _parseSingleEnd(plasmidFastqFile,spHandles,
                                              bcHandles,edit_dist)
    barcodeSpacerCounts = _keepTargetSpacers(barcodeSpacerCounts,spacer_file,
                                             edit_dist = 2)
    
    barcodeCounts = _sumDict(barcodeSpacerCounts)
    if edit_dist > 0:
        barcodeSpacerCounts,barcodeCounts = collapseBarcodesExact(barcodeSpacerCounts,
                                                                  edit_dist,barcodeCounts)
    
    barcodeSpacerCounts = _thresholdBarcodes(barcodeSpacerCounts,threshold)
    barcodeToSpacer = _callBarcodes(barcodeSpacerCounts,0.95,True)
    
    barcodeCounts = _sumDict(barcodeSpacerCounts) 
    usedReads = 0
    for v in barcodeCounts.values():
        usedReads += v

    print 'Used reads = ' + str(usedReads)
    print 'Number of barcodes = ' + str(len(barcodeSpacerCounts))
    
    guides = {} 
    with open(spacer_file) as f:
        for line in f:
            fields = str.strip(line).split('\t')
            guides[fields[1]] = fields[0]
    
    barcodeToGene = {k:guides[v].split('_')[0] for k,v in barcodeToSpacer.items()}
    
    return barcodeToGene


# Input: dictionary of counters
# Output: dictionary with counters summed
def _sumDict(dictOfCounters):
    totalCounts = Counter() 
    for k,v in dictOfCounters.items():
        totalCounts[k] += sum(v.values()) 
    return totalCounts


# Input: name of plasmid fastq file, sequence before spacer, sequence before barcode, edit distance
# Output: dictionary with barcode as key, spacer as value 
def _parseSingleEnd(plasmidFastqFile, spHandles, bcHandles, edit_dist):
    totalReads = 0
    barcodeSpacerCounts = defaultdict(Counter) # Number of reads per barcode

    f = gzip.open(plasmidFastqFile,'rb')
    for read in SeqIO.parse(f,'fastq'):
        totalReads += 1
        
        # Find start and end positions of spacer & barcode sequences
        spStart = list(approxHammingSearch(spHandles[0],read.seq))
        spEnd = list(approxHammingSearch(spHandles[1],read.seq))
        bcStart = list(approxHammingSearch(bcHandles[0],read.seq))
        #bcEnd = list(approxHammingSearch(bcHandles[1],read.seq))

        # Ensure that the spacer & barcode handles only map to one position on the read
        if len(spStart) >= 1 and len(spEnd) >= 1 and len(bcStart) >= 1:
            spStart = spStart[0][1] + len(spHandles[0]) 
            spEnd = spEnd[-1][1]
            bcStart = bcStart[0][1] + len(bcHandles[0])
            #bcEnd = bcEnd[-1][1]
            

            # Ensure the barcode is exactly 12 bp
            # if bcEnd-bcStart == 12:
            spacer = str(read.seq[spStart:spEnd])
            barcode = str(read.seq[bcStart:bcStart+BARCODE_LENGTH])
            
            quals = read.letter_annotations['phred_quality']
            spacerQuals = quals[spStart:spEnd]
            barcodeQuals = quals[bcStart:bcStart+BARCODE_LENGTH]

            if len(barcode) == BARCODE_LENGTH and belowQual(spacerQuals,20) <= 3 and \
                    belowQual(barcodeQuals,20) <= 3:
                barcodeSpacerCounts[barcode][spacer] += 1 
    f.close()
    return barcodeSpacerCounts
     

# Input: name of plasmid fastq file, sequence before spacer, sequence before barcode, edit distance
# Output: dictionary with barcode as key, spacer as value 
def _parsePairedEnd(plasmidFastqFile,plasmidFastqFileMate,spHandles,bcHandles,edit_dist):
    bcLen = len(bcHandles[0])
    spLen = len(spHandles[0])
    
    # Parse spacer fastq file (read 1)
    spacers = defaultdict(str)
    if plasmidFastqFile.split('.')[-1] == 'gz':
        f = gzip.open(plasmidFastqFile)
    else:
        f = open(plasmidFastqFile)
    for read in SeqIO.parse(f,'fastq'):
        spStart = list(approxHammingSearch(spHandles[0],read.seq))
        if len(spStart) >= 1:
            spStart = argmin(spStart) + spLen
            spacer =  str(read.seq[spStart:spStart+20])
            quals = read.letter_annotations['phred_quality']
            spQuals = quals[spStart:spStart + 20]
            if belowQual(spQuals,20) < 2:
                spacers[str(read.id)] = spacer 
    f.close()
    
    # Parse barcode fastq file (read 2)
    barcodeSpacerCounts = defaultdict(Counter)
    totalReads = 0
    
    if plasmidFastqFile.split('.')[-1] == 'gz':
        f = gzip.open(plasmidFastqFileMate)
    else:
        f = open(plasmidFastqFileMate)
    for read in SeqIO.parse(f,'fastq'):
        totalReads += 1
        readSeq = read.seq.reverse_complement()
        bcStart = list(approxHammingSearch(bcHandles[0],readSeq))
        bcEnd = list(approxHammingSearch(bcHandles[1],readSeq))
         
        if len(bcEnd) >= 1 and len(bcStart) >= 1:
            bcStart = argmin(bcStart) + len(bcHandles[0]) 
            bcEnd = argmin(bcEnd)
            barcode = str(readSeq[bcStart:bcEnd])
            quals = read.letter_annotations['phred_quality']
            barcodeQuals = quals[bcStart:bcEnd]
            if (len(barcode) == BARCODE_LENGTH or len(barcode) == BARCODE_LENGTH - 1) \
                    and belowQual(barcodeQuals,20) < 2 and str(read.id) in spacers:
                spacer = spacers[str(read.id)]
                barcodeSpacerCounts[barcode][spacer] += 1
    f.close()
    return barcodeSpacerCounts 


# Input: barcode-spacer mapping dictionary, minimum read threshold
# Output: barcode-spacer mappings with all mappings below minreads removed
def _thresholdBarcodes(barcodeSpacerCounts, minreads):
    for bc,spacerDict in barcodeSpacerCounts.items():
        for spacer,count in spacerDict.items():
            if count < minreads:
                del barcodeSpacerCounts[bc][spacer]
        
        if len(barcodeSpacerCounts[bc]) == 0:
            del barcodeSpacerCounts[bc]
    
    return barcodeSpacerCounts

# Input: barcode-spacer mapping counts, min read fraction
# Output: barcode-spacer mappings if one barcode-spacer mapping contains at
#         least the min read fraction of the reads among all possible mappings
def _callBarcodes(barcodeSpacerCounts,minFraction=0.9,print_collision=False):
    barcodeToSpacer = defaultdict(str)
    totalBarcodes = len(barcodeSpacerCounts)
    for barcode,counts in barcodeSpacerCounts.items():
        totalCounts = sum(counts.values())
        collision = True
        for key,count in counts.items():
            if float(count)/totalCounts >= minFraction:
                barcodeToSpacer[barcode] = key
                collision = False
                break
        
        if collision:
            del barcodeSpacerCounts[barcode]
        
    collisionRate = float(totalBarcodes-len(barcodeToSpacer))/totalBarcodes
    if print_collision:
        print 'Collision Rate = ' + str(collisionRate)

    return barcodeToSpacer


# Input: barcode-spacer mappings; file with all library spacer sequences, edit_dist
# Output: barcode-spacer mappings with only spacers that are in the library
def _keepTargetSpacers(barcodeSpacerCounts,libraryGuideFile,edit_dist=2):
    guides = set()
    with open(libraryGuideFile) as f:
        for line in f:
            guides.add(str.strip(line.split('\t')[1]))
    
    for barcode,spacerCounts in barcodeSpacerCounts.items():
        newSpacerCounts = Counter()
        for spacer in spacerCounts:
            guideSpacer = exactExists(guides,spacer,edit_dist)
            if guideSpacer:
                if type(guideSpacer) == str:
                    newSpacerCounts[guideSpacer] += barcodeSpacerCounts[barcode][spacer] 
                else:
                    newSpacerCounts[guideSpacer[0]] += barcodeSpacerCounts[barcode][spacer]

        if len(newSpacerCounts) == 0:
            del barcodeSpacerCounts[barcode]
        else:
            barcodeSpacerCounts[barcode] = newSpacerCounts

    return barcodeSpacerCounts


if __name__ == '__main__': main()
