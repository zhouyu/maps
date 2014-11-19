#!/usr/bin/env python
"""
Get intervals around unique TES (Transcription End Site)
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

    parser.add_option("--len_up", dest="len_up", default=1,
        help="length upstream of TES[default 1]", type="int")

    parser.add_option("--len_dn", dest="len_dn", default=0,
        help="length downstream of TES[default 0]", type="int")

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

    for gr in reader_gene(open(f_gene), 'bed'):
        if gr.get_strand() == "+":
            tes = gr.get_txEnd()
            s = tes - options.len_up
            e = tes + options.len_dn
        else:
            tes = gr.get_txStart()
            s = tes - options.len_dn
            e = tes + options.len_up

        outhandle.write("\t".join(map(str, [
            gr.get_chrom(), s, e, gr.get_score(), tes, gr.get_strand()]))+"\n")

    if options.outfile:
        outhandle.close()

if __name__ == '__main__':
    main()

