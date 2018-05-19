#!/bin/bash

cellBarcodePath=
fastqName=
outputGenotypeDict=
barcodeToGene=

# Convert from fastq to bam	
java -Xmx8g -jar $HOME/Yan_SSD/PicardTools/picard.jar FastqToSam F1=$fastqName\_R2.fastq.gz O=$fastqName\_R2.bam SM=$fastqName SORT_ORDER=unsorted TMP_DIR=tmp/

# Tag Read2 bam with Read1 cell/molecule barcodes
python ~/genotyping-matrices/convert_gRNA.py $fastqName\_R1.fastq.gz $fastqName\_R2.bam 
samtools index $fastqName\_R2.tagged.bam

# Parse cell barcodes and store in python dictionary
python $HOME/genotyping-matrices/parseCellBarcodes.py $cellBarcodePath $fastqName\_R2.tagged.bam 20 'GGCTGTTACGCG' 'CTACTGAC'
 
# Plot distribution of reads per UMI in order to set appropriate thresholds
python $HOME/genotyping-matrices/plotUMIdist.py $fastqName\_R2_cell_barcodes.pickle 0.01

# Genotype cells: outputs a genotype dictionary in csv format
python $HOME/genotyping-matrices/genotypeCells.py $fastqName\_R2_cell_barcodes.pickle $outputGenotypeDict $barcodeToGene 0.1 0.1 0.01
