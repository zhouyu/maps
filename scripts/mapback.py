#!/usr/bin/env python
"""
Mapping location from relative to absolute location

"""

import os
import sys
from maps.backmap import Backmapper


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Usage: %s infile.bed relative.bed outfile.bed' % sys.argv[0]
        sys.exit()
        
    f_in = sys.argv[1]
    f_relative = sys.argv[2]
    assert os.path.exists(f_in)
    assert os.path.exists(f_relative)
    outfile = sys.argv[3]
    
    mapper = Backmapper(f_relative)
    outhandle = open(outfile, 'w')
    mapper.backmap_file(f_in, outhandle=outhandle)
    outhandle.close()
    
    
