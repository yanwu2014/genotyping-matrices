import cPickle as cp
import numpy as np
from functools import partial
import time
import subprocess as sp
import pandas as pd
import sys
import shlex
from collections import defaultdict
import copy

# Represents a single exon
class Exon:

    def __init__(self,range,number,gene_id):
        self.gene_id = gene_id
        self.range = range
        self.length = np.abs(range[0]-range[1])
        self.number = number

# Represents a single gene
class Gene:

    def __init__(self,chr,range,id,name,attributes={},metadata={}):
        self.chr = chr # Chromosome
        self.name = name
        self.range = range # Gene start and ending positions
        self.id = id # Gene ID 
        self.exons = []
        self.attributes = attributes # Additional attributes in dict format
        self.metadata = metadata

# Load gencode gtf file as annotations
def loadAnnotationGTF(filename):
    features = {}
    
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            
            fields = str.strip(line).split('\t')
            start = int(fields[3])
            end = int(fields[4])
            chr = fields[0]
            
            metadata = {}
            for field in fields[8].split(';'):
                if str.strip(field):
                    key,val = shlex.split(field)
                    metadata[str.strip(key)] = str.strip(val)
            
            gene_id = str.strip(metadata['gene_id'])
            gene_id = gene_id.split('.')[0]
            gene_type = metadata.get('gene_type',None)
                
            gene_name = metadata.get('gene_name',None)  
            attr = {'name':gene_name, 'type':gene_type}

            if fields[2] == 'gene':
                feature = Gene(chr, (start,end), gene_id, gene_name,
                               attributes=attr, metadata=metadata)
                features[gene_id] = feature
            
            elif fields[2] == 'exon':
                exon_number = int(metadata['exon_number']) 
                
                exon = Exon((start,end),exon_number,gene_id)
                length = np.abs(start-end)
                if gene_id in features:
                    features[gene_id].exons.append(exon)
                    
                else:
                    raise ValueError('GTF format wrong: orphaned exon')
                 
    return features

