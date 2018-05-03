from datetime import datetime
from sys import argv
import copy
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os
import itertools
import statistics as s

def read_mapping_file(map_file):
    """TODO: Docstring for read_mapping_file.

    :map_file):: TODO
    :returns: dict with format of names 

    OLD - NEW
    """
    pairs = {}
    with open(map_file, 'r') as f:
        for line in f:
            clean = line.strip().split()
            pairs[clean[0]] = clean[1]
    return pairs

def GenomeReader(GenomeFile):
    """
    Arg: Takes in Genome File
    Rtrns: Returns a dictionary, Genome Scaffolds. 
    
    Keys genomic scaffold names being the
    keys - and the actual sequence being the value. 
    """
    GenomeScaffolds = {}
    with open(GenomeFile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                NamedSeq = line
                GenomeScaffolds[NamedSeq] = ""
            else:
                GenomeScaffolds[NamedSeq] += line
    return GenomeScaffolds


def rename_dict(mapping_dict, fasta_dict):
    """TODO: Docstring for rename_dict.

    :mapping_dict: TODO
    :fasta_dict: TODO
    :returns: TODO

    """

    for genename, seq in fasta_dict.items():
        isoalte_old_name = genename.split(' ')
        seq_name = isoalte_old_name[0].replace('>','')

        new_name = ">" + mapping_dict[seq_name]

        combine_new_string = new_name + ' ' + ' '.join(isoalte_old_name[1:])
        print(combine_new_string)
        print(seq)

        
read_in_genome = GenomeReader(argv[1])
map_dict = read_mapping_file(argv[2])
rename_dict(map_dict,read_in_genome)


def get_parser():
    parser = argparse.ArgumentParser(description='Software to read in fasta \
            file and do basic functionality that is often required  ')
    parser.add_argument('-f','--fasta', help='fasta file to read', \
            required=True, dest='f')
    parser.add_argument('-m','--map', help='reads in mapping file',\
            required=False,dest='s') 
    #Must add output file
    #parser.add_argument('-l','--locsplit', help='location in ccaffold to split',\
    #        required=False,dest='l') 


if __name__ == "__main__":
    args = get_parser().parse_args()
    StartTime = datetime.now()

    read_in_genome = GenomeReader(args.f)
    map_dict = read_mapping_file(args.m)
    rename_dict(map_dict,read_in_genome)
   
  
endtime = datetime.now()
finaltime = endtime - StartTime 

print ("Total Time %s" % (finaltime))

