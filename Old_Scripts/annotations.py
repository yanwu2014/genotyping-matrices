'''
Data structure & code to deal with gene and disease annotations
'''

from operator import itemgetter
import cPickle as cp
import sys

# Main can be used for testing.
def main():
    ccdsAnnotation = GeneAnnotation(sys.argv[1],'ccds')
    with open('ccds_annotations.pickle','w') as f:
        cp.dump(ccdsAnnotation,f)
    

# Gene object class represents a single gene
class Gene:
    
    def __init__(self,id):
        self.id = id # Gene id
        self.chr = None # Chromosome
        self.name = None # HUGO name
        self.codingStart = None # Coding start
        self.codingEnd = None # Coding end
        self.exons = [] 

# Object representing an group of genes
class GeneAnnotation:
    
    # Construct GeneAnnotation from either RefSeq or CCDS file
    def __init__(self, file, filetype):
                
        self.genes = {} # Genes represented as a dictionary indexed by gene name
             
        if str.lower(filetype) == 'ccds':
            self.reference = 'ccds'
            self._load_ccds_(file)
        
        elif str.lower(filetype) == 'refseq':
            self.reference = 'refseq'
            self._load_refseq_(file)
    
        else:
            raise ValueError('Invalid Filetype')
    
    # Load ccds annotation file
    def _load_ccds_(self, ccds_file):
        with open(ccds_file,'r') as f:
            f.readline().split('\t')
            for line in f:
                line = line.split('\t')
                status = str.strip(line[5])
                
                if status == 'Public':
                    chr = 'chr' + str.strip(line[0])
                    name = str.strip(line[2])
                    id = str.strip(line[4])
                    codingStart = int(line[7])
                    codingEnd = int(line[8])

                    exons = line[9].strip('[]').split(', ')
                    exons = [(int(x.split('-')[0]),int(x.split('-')[1])) \
                             for x in exons]
                    
                    thisGene = Gene(id)
                    thisGene.name = name
                    thisGene.chr = chr
                    thisGene.codingStart = codingStart
                    thisGene.codingEnd = codingEnd
                    thisGene.exons = exons
                    
                    self.genes[name] = thisGene
                        

    # Load refSeq annotations 
    def _load_refseq_(self, refSeq_file):
        with open(refSeq_file,'r') as f:
            f.readline().split('\t')
            for line in f:
                line = line.split('\t')
                
                id = str.strip(line[0])
                name = str.strip(line[1])
                chr = str.strip(line[2])
                
                exonStarts = str.strip(line[7]).split(',')
                exonEnds = str.strip(line[8]).split(',')
                del exonStarts[-1]
                del exonEnds[-1]
                exonStarts = [int(x) for x in exonStarts]
                exonEnds = [int(x) for x in exonEnds]
                exons = zip(exonStarts,exonEnds)
                
                codingStart = int(line[5])
                codingEnd = int(line[6])

                thisGene = Gene(id)
                thisGene.name = name
                thisGene.chr = chr
                thisGene.codingStart = codingStart
                thisGene.codingEnd = codingEnd
                thisGene.exons = exons

            self.genes[name] = thisGene


    def getExons(self,id):
        if id in self.genes:
            return self.genes[id].exons
        else:
            raise ValueError('ID not in reference database')
    
    def getChr(self,id):
        if id in self.genes:
            return self.genes[id].chr
        else:
            raise ValueError('ID not in reference database')

    def getName(self,id):
        if id in self.genes:
            return self.genes[id].name
        else:
            raise ValueError('ID not in reference database')
    
    

if __name__ == '__main__': main()
