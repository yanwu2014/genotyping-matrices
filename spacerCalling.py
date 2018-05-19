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
import spacerCalling as sc
import subprocess as sp

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


# Input: dataFrame, cellBarcodes object, barcode to gene mapping
# Output: genotyped dataframe
def callBarcodes(cellBarcodes, barcodeToGene, readThresh,
                 umiThresh, edit_dist = 1, minUMIFrac = 0.0, minUMIReads = 0, 
                 minGBCreads = 5, minGBCumis = 2):
    
    for cell_bc in cellBarcodes:
        if minUMIFrac > 0:
            cellBarcodes[cell_bc].countBarcodes(minFrac = minUMIFrac)
        else:
            cellBarcodes[cell_bc].countBarcodes(minReads = minUMIReads)
        
        cellBarcodes[cell_bc].callBarcode(barcodeToGene, readThresh, umiThresh,
                                          minGBCreads, minGBCumis)
    
    genotype_dict = defaultdict(list)
    for cell_bc,cell_obj in cellBarcodes.items():
        if cell_obj.type == 'usable':
            for genotype in cell_obj.genotype:
                genotype_dict[genotype].append(cell_bc)
    
    return genotype_dict

    


# Input: mapped & corrected BAM file from DropSeq, barcodeToSpacer mapping dict, 
#        cell readcount file from DropSeq, number of desired core barcodes
# Output: dict with cell barcodes as keys, Cell objects as vals 
def getCellBarcodes(bamFile, coreCells, barcode_length, bc_start_handle,
                    bc_end_handle, edit_dist = 1):
    
    bamFile = ps.AlignmentFile(bamFile, 'rb', check_sq=False) # load bam file
    
    cellBarcodes = {} 
    nreads = 0
    for read in bamFile.fetch(until_eof=True):
        molTag = _getTag(read,'XM')
        cellTag = _getTag(read,'XC')
        
        guideBarcode = _getBarcode(read, barcode_length, bc_start_handle, bc_end_handle)
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
def _getBarcode(read, barcode_length, bc_start_handle, bc_end_handle):
    # Make sure read maps to end of reference (right before barcode)
    #if len(read.query_sequence) < 65:
    #    return False
    
    # ensure that the sequence right before the barcode is what we expect
    bcStart = list(approxHammingSearch(bc_start_handle, read.query_sequence))
    if len(bcStart) < 1: 
        return False
    left_pointer = argmin(bcStart) + len(bc_start_handle)
    
    bcEnd = list(approxHammingSearch(bc_end_handle, read.query_sequence))
    if len(bcEnd) < 1: 
        return False
    right_pointer = argmin(bcEnd)

    barcode = read.query_sequence[left_pointer:right_pointer]

    # Ensure the read covers the entire barcode
    if (len(barcode) < barcode_length - 1) or (len(barcode) > barcode_length + 1):
        return False

    quals = read.query_qualities[left_pointer:right_pointer]
    if belowQual(quals, 20) > 1:
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


if __name__ == '__main__': main()
