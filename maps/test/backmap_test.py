#!/usr/bin/env python

import sys
import unittest
from StringIO import StringIO
from rp3seq.backmap import Ival, GeneRna, Backmapper

###############################################################################

class Backmapper_Test(unittest.TestCase):
    """Test """
    def setUp(self):
        self.f_gene = 'data/backmap_gene.bed'
        self.f_read = 'data/backmap_tag.bed'
        self.mapper = Backmapper(self.f_gene)
        
    def test_backmap(self):
        out = StringIO()
        self.mapper.backmap_file(self.f_read, out)
        expected = [
['chr1', '380', '400', '1', '0', '+', '380', '380', '0,0,0', '1', '20', '0'],
['chr2', '120', '140', '2', '2', '-', '120', '120', '0,0,0', '1', '20', '0'],
['chr1', '190', '310', '3', '3', '+', '190', '190', '0,0,0', '3', '10,10,10', '0,20,110'],
['chr2', '300', '330', '4', '4', '-', '300', '300', '0,0,0', '1', '30', '0']]
        obs = [line.split("\t") for line in out.getvalue().rstrip().split("\n")]
        assert expected == obs

############################################################################### 

class GeneRna_Test(unittest.TestCase):
    def setUp(self):
        self.e1 = Ival(100, 200)
        self.e2 = Ival(210, 220)
        self.e3 = Ival(300, 420)
        self.grna1 = GeneRna('chr1', '+', (self.e1, self.e2, self.e3))
        self.grna2 = GeneRna('chr2', '-', (self.e1, self.e2, self.e3))
        
    def test__create_blocklist(self):
        assert self.grna1._create_blocklist() == [Ival(0, 100), Ival(100, 110), Ival(110, 230)]
        assert self.grna2._create_blocklist() == [Ival(0, 120), Ival(120, 130), Ival(130, 230)]
        
    def test_backmap(self):
        assert self.grna1.backmap_ival(190, 210) == [Ival(380, 400)]
        assert self.grna2.backmap_ival(190, 210) == [Ival(120, 140)]
        assert self.grna1.backmap_ival(90, 120) == [Ival(190, 200), Ival(210, 220), Ival(300, 310)]
        assert self.grna2.backmap_ival(90, 120) == [Ival(300, 330)]
    
    def test_as_bed12(self):
        assert self.grna1.as_bed12('uc1') == (
            'chr1', 100, 420, 'uc1', 0, '+', 100, 100, '0,0,0', 3, '100,10,120', '0,110,200')
        assert self.grna2.as_bed12('uc2') == (
            'chr2', 100, 420, 'uc2', 0, '-', 100, 100, '0,0,0', 3, '100,10,120', '0,110,200')
        #print "\t".join(map(str, self.grna1.as_bed12('uc1')))
        #print "\t".join(map(str, self.grna2.as_bed12('uc2')))
        
###############################################################################

if __name__ == '__main__':
    unittest.main()
    
    
