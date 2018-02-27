import sys
import cPickle as cp
from spacerCalling import *
import pandas as pd

cell_barcodes_file = sys.argv[1]
barcode_gene_file = sys.argv[2]

with open(cell_barcodes_file) as f:
    cell_barcodes = cp.load(f)

barcodeToGene = {}
with open(barcode_gene_file) as f:
    for line in f:
        fields = str.strip(line).split("\t")
        barcodeToGene[fields[0]] = fields[1]

cell_umi_counts_dict = {}
for cbc,cell in cell_barcodes.items():
    cell.countBarcodes()
    guideUMICounts = Counter()
    for gbc in cell.gbcUMICounts:
        libGBC = exactExists(barcodeToGene, gbc, edit_dist = 1)
        if libGBC and type(libGBC) == str:
            guide = barcodeToGene[libGBC]
            guideUMICounts[guide] += cell.gbcUMICounts[gbc]

    cell_umi_counts_dict[cbc] = guideUMICounts
    
cell_umi_counts = pd.DataFrame.from_dict(cell_umi_counts_dict)
cell_umi_counts = cell_umi_counts.fillna(0.0)

out_file = cell_barcodes_file.replace(".pickle", ".tsv")
cell_umi_counts.to_csv(out_file, sep = "\t")
