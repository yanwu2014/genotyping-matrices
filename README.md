# genotyping-matrices

Scripts and python modules for genotyping cells from perturbation screens followed by scRNA-seq. Currently compatible with the 10X genomics Chromium V2 chemistry.

Required input: 
* Genotype barcodes amplified from 10X genomics scRNA-seq cDNA: Read1 will be the cell/molecule barcode and Read2 should have the genotype barcode. Read2 needs to be in bam format. You can use PicardTools FastqToBam to convert the Fastq to bam, or align Read2 to the known genotype scaffold as an initial filtering step.
* Text file mapping barcodes to genotypes. These genotypes can be genes, CRISPR gRNAs, or any other genotype (mutations, etc.)

Usage (Replace filenames in brackets with the names of your files):
1. Tag Read2 genotype reads with cell and molecule barcodes using convert_gRNA.py: `python convert_gRNA.py [read1_cell_molecule_barcodes].fastq [read2_genotype_barcodes].bam`
2. Parse the cell barcodes with parseCellBarcodes.py. You'll need to edit the parseCelLBarcodes.py script and replace BC_START_HANDLE and BC_END_HANDLE with the sequences upstream and downstream of your genotype barcode, and BARCODE_LENGTH with the length of your genotype barcode. `python parseCellBarcodes.py `
