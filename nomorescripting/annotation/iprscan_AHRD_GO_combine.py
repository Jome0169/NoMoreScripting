#!/usr/local/bin/python3
import getopt
import sys
import os
import re
import copy
from sys import argv
from datetime import datetime
import argparse

def read_ahrd_csv(ahrd_file):
    """reads in ahrd file. Parses both the gene name and the funciton

    :ahrd_file: TODO
    :returns: TODO

    """
    annotation_key = {}

    with open(ahrd_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split('\t')
            if "#" in clean_line[0]:
                pass
            elif clean_line == None:
                pass
            else:
                annotation_key[clean_line[0]] = clean_line[1:]
    return annotation_key

def read_iprscan(iprscan_file):
    """reads in interproscan output file. Creates a list for later processsing.

    :iprscan_file: interpro scan file to read
    :returns: a nested list.

    """
    ipscan_dict = {}
    with open(iprscan_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split()
            if clean_line[0] not in ipscan_dict:
                ipscan_dict[clean_line[0]] = [clean_line]
            elif clean_line[0] in ipscan_dict:
                ipscan_dict[clean_line[0]].append(clean_line)
    return ipscan_dict


def read_GO_file(go_file):
    """reads in GO file, returns a list like element

    :go_file: TODO
    :returns: TODO

    """
    go_list = []

    with open(go_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split()
            go_list.append(clean_line)
    return go_list

def read_fasta_file(fasta_file):
    """Read in headers from fasta file. Creates a dictioninary for genes so
    that functional inforamtion can be added to the value inforamtion

    :fasta_file: fasta file
    :returns: dictionary with genes names as key and empty list as value

    """
    header_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                clean_header = line.strip().replace('>','')
                header_dict[clean_header] = []
            else:
                pass
    return header_dict

def combine_dict_AHRD(key_dict, ahrd_lines):
    """begins by appending ahrd information to the appropriate gene value in
    the gene dict

    :key_dict: TODO
    :ahrd_lines: TODO
    :returns: TODO

    """
    for item in ahrd_lines:
        if ahrd_lines[0] in key_dict:
            create_ahrd_list = [ahrd_lines[1], ahrd_lines[3]]
            key_dict[ahrd_lines[0]].append(ahrd_lines[1])
            key_dict[ahrd_lines[0]].append(ahrd_lines[1])
    return key_dict

def combine_dict_iprscan(key_dict,ipscan_dict):
    """Appends ipscan information to the gene dict

    :key_dict: TODO
    :ipscan_dict: TODO
    :returns: TODO

    """
    for ipscan_key, ipscan_info in ipscan_dict.items():
        if ipscan_key in key_dict:

            key_dict[ipscan_key].append(ipscan_info)
    return key_dict


def digest_iprscan_values(key_dict):
    """TODO: Docstring for digest_iprscan_values.

    :key_dict: TODO
    :returns: TODO

    """
    for gene_key, gene_val in key_dict.items():






def get_parser():
    parser = argparse.ArgumentParser(description='The purpose of this script it\
            to take in AHRD file files, IPRscan, as well as GO output and\
            combines them to create the final functional annoation xlxs file.')
    parser.add_argument('-a','--ahrd', help='AHRD output file', required=True, dest='a')
    parser.add_argument('-f','--fasta', help='fasta file', required=true, dest='p')
    parser.add_argument('-ip','--interpro', help='interproscan file to add to output', required=true, dest='ip')
    parser.add_argument('-go','--goterms', help='go term output file to add to AHRD output file', required=true, dest='go')
    parser.add_argument('-o','--output', help='File to Write output to. Will \
            Remove file if finds a conflicting name', required=True, dest='o')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    StartTime = datetime.now()
    args = get_parser().parse_args()





#Speed Things
EndTime = datetime.now()
FinalTime = EndTime - StartTime

print ("Total Time %s" % (FinalTime))







