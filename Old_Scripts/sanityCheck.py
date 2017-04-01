# Check to make sure cDNA barcodes correspond to gRNA barcodes

import spacerCalling as sc
import cPickle as cp
import sys

def main():
    checkBarcodes()

def checkBarcodes():
    gRNAbamFile = sys.argv[1]
    cDNAbamFile = sys.argv[2]
    cellReadCounts = sys.argv[3]
    barcodeToGeneFile = sys.argv[4]
    numCells = int(sys.argv[5])

    with open(barcodeToGeneFile) as f:
        barcodeToGene = cp.load(f)
    
    gRNAcellBarcodes = sc.getCellBarcodes(gRNAbamFile,cellReadCounts,numCells,barcodeToGene)
    cDNAcellBarcodes = sc.getCellBarcodes(cDNAbamFile,cellReadCounts,numCells,barcodeToGene)
    
    for barcode,cDNA_cell in cDNAcellBarcodes.items():
        if barcode in gRNAcellBarcodes:
            gRNA_cell = gRNAcellBarcodes[barcode]
            print 'gRNA barcodes: '
            for barcode,count in gRNA_cell.barcodeReadCounts.items():
                print barcode + '\t' + str(count)
            
            print 'cDNA barcodes: '
            for barcode,count in cDNA_cell.barcodeReadCounts.items():
                print barcode + '\t' + str(count)

            print '--------------------------------------------\n'

if __name__ == '__main__': main()
