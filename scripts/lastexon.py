#!/usr/bin/env python
"""
Get length of unique last exons
"""

import os
import sys
import optparse
from maps.io_utils.gene import reader_gene

def process_command_line(argv):
    if argv is None:
        argv = sys.argv[1:]
    usage = "usage: %prog [options] gene.bed"
    parser = optparse.OptionParser(usage,
        formatter=optparse.TitledHelpFormatter(width=178),
        add_help_option=None)

    parser.add_option('-q', '--quiet',
        action='store_false', dest='verbose', default=True,
        help="don't print status messages to stdout")

    parser.add_option("-o", "--outfile", dest="outfile",
        help="out file name", type="str")

    parser.add_option('-h', '--help', action='help',
        help='Show this help message and exit.')

    (options, args) = parser.parse_args(argv)
    if len(args) < 1:
        parser.error('No required parameters')

    return options, args

def main(argv=None):
    options, args = process_command_line(argv)
    f_gene = args[0]
    assert os.path.exists(f_gene)
    outhandle = sys.stdout
    if options.outfile:
        outfile = options.outfile
        outhandle = open(outfile, "w")

    llexons = set()
    for gr in reader_gene(open(f_gene), 'bed'):
        chrom = gr.get_chrom()
        exons = gr.get_exons()
        strand = gr.get_strand()
        if strand == "+":
            llexons.add((chrom, exons[-1][0], exons[-1][1], strand))
        else:
            llexons.add((chrom, exons[0][0], exons[0][1], strand))

    i = 0
    for chrom, s, e, strand in sorted(llexons):
        i += 1
        outhandle.write("\t".join(map(str, [
            chrom, s, e, i, e - s, strand]))+"\n")

    if options.outfile:
        outhandle.close()

if __name__ == '__main__':
    main()

