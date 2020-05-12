#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 08 May 2018 15:14:43

import sys, os, argparse, string
from Bio import SeqIO

from collections import Counter

def ParseFastq():
    ''' main scripts '''
    ### Move the first 3 letters into the fastq header
    for read in SeqIO.parse(sys.argv[1], 'fastq'):
        print "@" + read.id + "|"+ read.seq[0:3]
        print read.seq[4:]
        print "+"
        print "".join([chr(x+33) for x in read.letter_annotations['phred_quality']])[4:]


def main():
    '''main scripts
    '''
    ParseFastq()


if __name__ == "__main__":
    main()

