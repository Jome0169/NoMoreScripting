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
    parser.add_argument('-s','--scaffold', help='scaffold to split',\
            required=False,dest='s') 
    parser.add_argument('-l','--locsplit', help='location in ccaffold to split',\
            required=False,dest='l') 


if __name__ == "__main__":
    args = get_parser().parse_args()
    StartTime = datetime.now()

    file_dict = read_in_fasta(args.f)
    
    ##At some point these need to be fixed
    if args.s != None and args.l == None:
        #exit if split given without where
        print("If splitting file need location so split scaffold")
        sys.exit(-2)
    elif args.s == None and args.l != None:
        print("If splitting file need scaffold to split")
        sys.exit(-2)

    elif args.ra != None and args.s != None and args.l != None:
        broken_scaffold = split_scaffold(file_dict, args.s, args.l)
        brokenscafwriter(args.o,file_dict,args.s,broken_scaffold,args.ra)
    elif args.ra != None and args.s == None and args.l == None:
        rename_scaf_writer(args.o,file_dict,args.ra)
    elif args.ra == None and args.s != None and args.l != None:
        broken_scaffold = split_scaffold(file_dict, args.s, args.l)
        brokenscafwriter(args.o,file_dict,args.s,broken_scaffold,None)

    #If we just want to unwrap fasta
    elif args.ra == None and args.s == None and args.l == None \
    and args.uw != None:

        flaten_fasta(file_dict,args.o)


    elif args.ra == None and args.s == None and args.l == None \
    and args.uw == None and args.fl1 != None and args.fl2 != None:
    
        final_list = filter_fasta_dict(file_dict,args.fl1, args.fl2)
        print(final_list)
        flaten_fasta(final_list, args.o)
    
    elif args.ra == None and args.s == None and args.l == None \
    and args.uw == None and args.fl1 == None and args.fl2 == None and args.plt \
    == None and args.st \
    != None:
        report_fasta_stats(file_dict)
    
    
    elif args.plt != None:
        before_imoprt = datetime.now()
        from plt_scaf import *
        after_import = datetime.now()
        
        print(str(after_import - before_imoprt))

        chrom_len_dict = calculate_chrom_len(file_dict)
        ns_window_dict = calc_ns_in_window(file_dict, 100000)
        hist_plot_ns(chrom_len_dict,ns_window_dict, 100000 )
    else:
        print("NO GO")
   
endtime = datetime.now()
finaltime = endtime - StartTime 

print ("Total Time %s" % (finaltime))






