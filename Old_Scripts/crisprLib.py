from annotations import *
import subprocess as sp
import cPickle as cp

# Represents a single spacer sequence
class Spacer:
    def __init__(self,seq):
        self.seq = seq
        self.pos = None
        self.score = None

# Class that represents a single gene in a library
class Gene:

    def __init__(self,name):
        self.name = name
        self.chr = None
        self.id = None
        self.spacers = []
    
    def addSpacer(self,seq,score):
        thisSpacer = Spacer(seq)
        thisSpacer.score = score
        i = 0
        if len(self.spacers) == 0:
            self.spacers.append(thisSpacer)
        else:
            while score < self.spacers[i].score:
                i += 1
                if i == len(self.spacers):
                    self.spacers.append(thisSpacer)
                    break
            self.spacers.insert(i,thisSpacer)

    def deleteSpacer(self,spacer):
        self.spacers.remove(spacer)


# Class that represents a CRISPR library
class Library:
    
    # Create CRISPR library from a list of genes
    def __init__(self,geneList,refFile,barcodes,scaffold):
        self.barcodes = barcodes
        self.scaffold = scaffold
        self.genes = {}
        geneList = set(geneList)
        added = set()
        with open(refFile,'r') as f:
            for line in f:
                line = line.split('\t')
                gene = str.strip(line[0])
                seq = str.strip(line[1])
                if len(line) >= 3:
                    score = float(str.strip(line[2]))
                else:
                    score = 0.0
                if gene in geneList:
                    if not gene in self.genes:
                        thisGene = Gene(gene)
                        thisGene.addSpacer(seq,score)
                        added.add(gene)
                        self.genes[gene] = thisGene
                    else:
                        self.genes[gene].addSpacer(seq,score)


        if geneList != added:
            print 'Genes not found:'
            for gene in (geneList - added):
                print gene
    

    # Checks spacer sequences to ensure no-multimapping and that all spacers
    # map to exons of the genes they belong to
    def checkPositions(self,genomeFile,annotation,p):
        with open('hg38.pickle','r') as f:
            genome = cp.load(f)
        
        # Print spacers to fasta file
        with open('temp_fasta_file.fa','w') as f:
            for gene in self.genes.values():
                for i,spacer in enumerate(gene.spacers):
                    f.write('>' + gene.name + ',' + str(i) + '\n')
                    f.write(spacer.seq + '\n')
        
        bwtCmd = 'bowtie -k 8 -v 0 --best --suppress 2,5,6,7,8 -p ' + \
                  str(p) + ' ' + genomeFile + ' --quiet -f ' + 'temp_fasta_file.fa'
        output = sp.check_output(bwtCmd, shell=True)
        output = output.split('\n')
        
        unmapped = {}
        mapped = 0
        for line in output:
            if str.strip(line):
                line = line.split('\t')
                geneName = str.strip(line[0]).split(',')[0]
                rank = int(str.strip(line[0]).split(',')[1])
                chr = str.strip(line[1])
                pos = int(str.strip(line[2]))
                
                exons = annotation.genes[geneName].exons
                inExon = False
                
                assert len(genome[chr][pos:pos+23]) == 23
                if genome[chr][pos+21:pos+23] == 'GG' or genome[chr][pos+21:pos+23] == 'CC':
                    for exonStart,exonEnd in exons:
                        if pos <= int(exonEnd) and pos >= int(exonStart):
                            inExon = True
                    if not inExon:
                        if geneName in unmapped:
                            unmapped[geneName].add(int(rank))
                        else:
                            unmapped[geneName] = set([int(rank)])
        
        mismapped_items = 0
        for indices in unmapped.values():
            mismapped_items += len(indices)
        print 'Number of mismapped guides: ' + str(mismapped_items)
        #for k,v in unmapped.items():
        #    print k + ': ' + str(v)

        for name,indices in unmapped.items():
            spacers = self.genes[name].spacers
            newSpacers = []
            for i in range(0,len(spacers)):
                if not i in indices:
                    newSpacers.append(spacers[i])
            if len(newSpacers) == 0:
                print name
                del self.genes[name]
            else:
                self.genes[name].spacers = newSpacers
        
        sp.call(['rm','temp_fasta_file.fa'])
    

    # Print spacer sequences in tab delimited format to outputFile
    def printSpacers(self,outputFile,libName):
        libSize = 0
        with open(outputFile,'w') as f:
            for gene in self.genes.values():
                for i,spacer in enumerate(gene.spacers):
                    outSeq = self.barcodes[0] + self.scaffold[0] + spacer.seq + \
                                self.scaffold[1] + reverseComplement(self.barcodes[1])
                    f.write(libName + '-' + gene.name + '-' + str(i) + '\t' + outSeq + '\n')
                    libSize += 1
        print 'Library Size = ' + str(libSize)
    

    # Keep at most the top k spacer sequences for each gene
    def pruneSpacers(self,k):
        for gene in self.genes.values():
            if k < len(gene.spacers):
                gene.spacers = gene.spacers[0:k]

    def addGene(self,gene):
        self.guides[gene.name] = gene

    def deleteGene(self,gene):
        if gene.name in self.genes:
            del self.genes[gene.name]

    def hasGene(self,guide):
        if gene.name in self.gene:
            return True
        else:
            return False

# Reverse complements string text
def reverseComplement(text):
    rev_comp = []
    for i in range(0, len(text)):
        if text[i] == 'A': rev_comp.append('T')
        elif text[i] == 'C': rev_comp.append('G')
        elif text[i] == 'T': rev_comp.append('A')
        elif text[i] == 'G': rev_comp.append('C')
        else: 
            print text[i]
            raise ValueError('Invalid Nucleotide')
    return ''.join(reversed(rev_comp))

