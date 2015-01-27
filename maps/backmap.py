#!/usr/bin/env python
"""
Module for back map coordinates in mRNA to the genome
"""
import sys
from maps.io_utils.bed import BedReader
from operator import itemgetter


def as_bed12(chrom, exons, strand, name, cdsStart=None, cdsEnd=None, score=0):
    """As bed12 formatted string"""
    exons = sorted(exons, key=itemgetter(0, 1))
    exonCount = len(exons)
    assert exonCount > 0
    txStart = exons[0][0]
    txEnd = exons[-1][1]
    exon_starts = [e[0] for e in exons]
    exon_ends = [e[1] for e in exons]
    block_sizes = [exon_ends[i] - exon_starts[i] for i in xrange(exonCount)]
    block_starts = [exon_starts[i] - txStart for i in xrange(exonCount)]
    if not cdsStart or not cdsEnd:
        cdsStart = txStart
        cdsEnd = txStart
        
    return (chrom, txStart, txEnd, name, score, strand,
        cdsStart, cdsEnd, '0,0,0', exonCount,
        ",".join(map(str, block_sizes)),
        ",".join(map(str, block_starts)),
        )


class Backmapper(object):
    """Map coordinates in mRNA to the genome"""
    def __init__(self, f_gene):
        self.f_gene = f_gene
        self.gindex = self._create_geneindex()
        
    def backmap_file(self, f_tag, outhandle=None):
        """Back map reads in one file"""
        if not outhandle:
            outhandle = sys.stdout
        for b in BedReader(open(f_tag)):
            gname = b.get_chrom()
            assert gname in self.gindex
            if b.get_strand() == '-':  # relative position only in fwd strand
                continue
            g = self.gindex[gname]
            ivals = [(ival.get_start(), ival.get_end()) for ival in 
                g.backmap_ival(b.get_start(), b.get_end())]
            
            outhandle.write("\t".join(map(str, 
                as_bed12(g.get_chrom(), ivals, g.get_strand(), b.get_name(), 
                    score=b.get_score())
                ))+"\n")
            
    def _create_geneindex(self):
        """Create geneid:GeneRna object index"""
        ginx = dict()
        for g in BedReader(open(self.f_gene)):
            chrom = g.get_chrom()
            strand = g.get_strand()
            start = g.get_start()
            end = g.get_end()
            
            # calculate exonic Ivals
            if g.get_nfields() <= 6:
                exons = [Ival(start, end)]
            else:
                exon_starts = [int(i) + start 
                    for i in g.get_blockStarts().strip(',').split(',')]
                blockSizes = [int(i) 
                    for i in g.get_blockSizes().strip(',').split(',')]
                exon_ends = [base + offset 
                    for base, offset in zip(exon_starts, blockSizes)]
                exons = [Ival(exon_starts[k], exon_ends[k]) 
                    for k in xrange(len(exon_starts))]
                
            ginx[g.get_name()] = GeneRna(chrom, strand, tuple(exons))
        
        return ginx
    

class GeneRna(object):
    """Gene with its spliced RNA"""
    def __init__(self, chrom, strand, exon_list):
        self.chrom = chrom
        self.strand = strand
        self.exon_list = exon_list
        self.blocklist = self._create_blocklist()
        self.numblock = len(self.blocklist)
        
    def get_chrom(self):
        return self.chrom
    
    def get_strand(self):
        return self.strand
        
    def backmap_ival(self, start, end):
        """Back map one interval"""
        k = self.find_firstblock_contains(start)
        
        newivals = list()
        for i in xrange(k, self.numblock):
            gs = self.blocklist[i].get_genomic_start()
            bs = self.blocklist[i].get_start()
            be = self.blocklist[i].get_end()
            
            # new mRNA-like position in current block
            qs = max(bs, start) 
            qe = min(be, end)
            
            # relative position in block + block start
            hs = qs - bs + gs
            he = qe - bs + gs
            if self.strand == '-':
                hs = be - qe + gs
                he = be - qs + gs
            
            if hs == he: # length 0
                continue
            
            newivals.append(Ival(hs, he))
            if end <= be:
                break
            
        if self.strand == '-':
            newivals.reverse()
        
        return newivals
    
    def find_firstblock_contains(self, pos):
        """Find first block contain given position"""
        i = self.numblock - 1
        while i >= 0:
            if self.blocklist[i].contains_pos(pos):
                break
            else:
                i -= 1
        return i
    
    def _create_blocklist(self):
        """Create mRNA block list"""
        blocks = list()
        exons = sorted(self.exon_list, key=lambda exon: exon.get_start())
        if self.strand == '-':
            exons.reverse()
        
        curr_mrnalen = 0
        for exon in exons:
            gs = exon.get_start()
            bs = curr_mrnalen
            be = curr_mrnalen + exon.get_len()
            blocks.append(Block(gs, bs, be))
            
            curr_mrnalen += exon.get_len()
        
        return blocks
    
    def as_bed12(self, name=None):
        return as_bed12(self.chrom, [
            (e.get_start(), e.get_end()) for e in self.exon_list],
            self.strand, name)
        

class Ival(object):
    """
    Basic interval with start and end
    """
    def __init__(self, start, end):
        assert start <= end, "start must be less than end"
        self.start = start
        self.end = end
        
    def __cmp__(self, other):
        return cmp(self.start, other.start) or cmp(self.end, other.end)
    
    def __hash__(self):
        return hash((self.start, self.end))
    
    def __repr__(self):
        return "Ival(%d, %d)" % (self.start, self.end)

    def __str__(self):
        return self.__repr__()

    def get_start(self):
        """start position"""
        return self.start
    
    def get_end(self):
        """end position"""
        return self.end
    
    def get_len(self):
        """length of the interval"""
        return self.end - self.start

    def overlaps(self, other):
        """whether overlaps with other interval"""
        if self.start >= other.end or self.end <= other.start:
            return False
        else:
            return True
    
    def write_pos(self):
        """return position string"""
        return '%d-%d' % (self.start, self.end)
    
    def contains_pos(self, pos):
        return self.start <= pos and self.end > pos
    

class Block(Ival):
    """A spliced block in one gene"""
    def __init__(self, gstart, bstart, bend):
        super(Block, self).__init__(bstart, bend) # (start, end) in spliced gene
        self.gs = gstart # genomic start
    
    def __str__(self):
        return super(Block, self).__str__() + "%d" % self.gs
    
    def get_genomic_start(self):
        """Return genomic start position"""
        return self.gs
        
