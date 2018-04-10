#!/usr/local/bin/python3
import getopt
import sys
import os
import re
import copy
from sys import argv
from datetime import datetime
import argparse




def read_in_gff(gff_file):
    """TODO: Docstring for read_in_gff.

    :gff_file: TODO
    :returns: TODO

    """
    all_gff_records = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            cleanline = line.strip().split()
            all_gff_records.append(cleanline)
    return all_gff_records



def only_print_genes(gff_list):
    """TODO: Docstring for only_print_genes.

    :gff_list: TODO
    :returns: TODO

    """
    for item in gff_list:
        if len(item) > 2:
            if item[2] == 'gene' or item[2] == 'mRNA' or item[2] == 'CDS' \
                or item[2] == 'exon' or item[2] == 'three_prime_UTR' \
                or item[2] == "five_prime_UTR":

                    print('\t'.join(item))
            else:
                pass
        elif len(item) == 2:
            print(' '.join(item))


def get_parser():
    parser = argparse.ArgumentParser(description='The purpose of this script is \
            print only gene encoding regions into a gff3 file instead of all\
            extra material as output by MAKER.')
    parser.add_argument('-g','--gff', help='gff file to extract genes back into', required=True, dest='g')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    StartTime = datetime.now()
    args = get_parser().parse_args()
    
    read_gff = read_in_gff(args.g)
    
    only_print_genes(read_gff)

#Speed Things
EndTime = datetime.now()
FinalTime = EndTime - StartTime




