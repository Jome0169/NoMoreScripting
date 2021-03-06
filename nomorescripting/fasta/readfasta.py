from datetime import datetime
from sys import argv
import copy
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys
import os
import itertools
import statistics as s


def read_in_fasta(arg1):
    """Reads in Fasta file using biopyton

    :arg1: TODO
    :returns: TODO

    """
    file_dict = SeqIO.index(arg1, "fasta")
    return file_dict

def split_scaffold(fasta_file,scaf_name,loc_split):
    """Takes in a fasta file and split it at the location given by loc_split
    variable which should be a number

    :fasta_file: TODO
    :scaf_name: TODO
    :loc_split: TODO
    :returns: TODO

    """
    finallist = []
    try:
        int(loc_split)
    except ValueError:
        print("%s not a int") % loc_split
        sys.exit(-2)

    if scaf_name in fasta_file:
        total_fasta_len = len(fasta_file[scaf_name])
        firsthalf = fasta_file[scaf_name][0:int(loc_split)]
        secondhalf = fasta_file[scaf_name][int(loc_split)+1:]
        
        fastastringfirst = str(scaf_name) + '_1'
        fastastringsecond = str(scaf_name) + '_2'

        finallist.append([fastastringfirst, firsthalf])
        finallist.append([fastastringsecond,secondhalf])
    
    elif scaf_name not in fasta_file:
        print("scaffold %s not found in file" % scaf_name)
        exit(-2)

    return finallist


def filter_fasta_dict(fasta_dict, filter1, filter2):
    """TODO: Docstring for filter_fasta_dict.

    :fasta_dict: TODO
    :filter1: TODO
    :filter2: TODO
    :returns: TODO

    """
    final_file_list = {}

    for scaf,record in fasta_dict.items():
        if len(record.seq) >= int(filter1) and len(record.seq) <= int(filter2):
            final_file_list[scaf] = record
        elif len(record.seq) <= int(filter1) or len(record.seq) >= int(filter2):
            pass
    return final_file_list


def report_fasta_stats(fasta_dict):
    """takes in the fasta sequence object and reports back some basics
    statistics. 

    Mean, Medien, SmallestLegnth, LongestLength,

    :fasta_dict: TODO
    :returns: TODO

    """
    total_val = 0
    number_over_25kb = 0
    store_lens = []
    for scaf_name, val in fasta_dict.items():
        take_seq_len = len(val.seq)
        total_val += take_seq_len
        store_lens.append(int(take_seq_len))
        if take_seq_len >= 25000:
            number_over_25kb += 1 
    store_lens.sort()
    mean_calc = s.mean(store_lens)
    medien_calc = s.median(store_lens)
    smalles_len = store_lens[0]
    largest_len = store_lens[-1]

    total_sum = sum(store_lens)/2
    
    
    #Print calculated stats
    print("The Smallest len is %s" % str(smalles_len))
    print("The largest len is %s" % str(largest_len))
    print("The medien len is %s" % str(medien_calc))
    print("The mean len is %s" % str(mean_calc))
    print("The total number of basepairs is %s" %str(total_val))

    print("The number of scaffolds over 25kB is %s" % str(number_over_25kb))
    print()


def brokenscafwriter(outputname,fasta_file,scaf_name,borken_scaf,scaf_rename):
    """Writes the broken scaffold to 2 outputs. Borken scaffold, rest of
    scaffolds, 

    :outputname: basename of the broken file output
    :fasta_file: input fasta dict from biopython
    :scaf_name: scaf name that was supposed to be broken
    :borken_scaf: list of broken scaf seq

    """
    #Broken scaffold file name
    broken_scaf_file_name = str(scaf_name) +  ".broken.fasta"
    
    #If renaming string, create variable indicating if string is real
    append_to_scaf_string = None
    if scaf_name == False:
        pass
    else:
        append_to_scaf_string = scaf_rename 

    #Remove Files if they exists
    try:
        os.remove(broken_scaf_file_name)
        os.remove(outputname)
    except OSError:
        pass
    
    #At some poitn probably make a function to go through and rename
    #dict keys based off stuff. To be done another time
    with open(outputname, 'a+') as f:
        for key, val in fasta_file.items():
            if key != scaf_name and append_to_scaf_string != None:
                f.write(">" + append_to_scaf_string + str(key))
                f.write('\n')
                f.write(str(val.seq))
                f.write('\n')
            elif key != scaf_name and append_to_scaf_string == None:
                f.write(">" + str(key))
                f.write('\n')
                f.write(str(val.seq))
                f.write('\n')
            else:
                pass
    
    with open(broken_scaf_file_name, 'a+') as z:
        for item in borken_scaf:
            if append_to_scaf_string != None:
                z.write(">" + append_to_scaf_string + item[0])
                z.write('\n')
                z.write(str(item[1].seq))
                z.write('\n')
            elif append_to_scaf_string == None:
                z.write(">" + item[0])
                z.write('\n')
                z.write(str(item[1].seq))
                z.write('\n')


def replace_name(fasta_file_dict, replace_file):
    """TODO: Docstring for replace_name.

    :fasta_file_dict: TODO
    :replace_file: TODO
    :returns: TODO

    """
    pairs = []
    with open(replace_file, 'r') as f:
        for line in f:
            cleanpair = line.split()
            if cleanpair == None:
                pass
            else:
                pairs.append(cleanpair)
    return pairs


def rename_scaf_writer(outputname,fasta_file,scaf_rename):
    """If appending string to scaffold name, write new output fasta file. This
    seemed easier than making a bunch of if then statments below.

    :outputname: TODO
    :fasta_file: TODO
    :scaf_rename: TODO
    :returns: TODO

    """
    #If renaming string, create variable indicating if string is real
    append_to_scaf_string = None
    if scaf_name == False:                                            
        pass
    else:
        append_to_scaf_string = scaf_rename 

    #Remove Files if they exists
    try:
        os.remove(outputname)
    except OSError:
        pass

    with open(outputname, 'a+') as f:
        for key, val in fasta_file.items():
            if append_to_scaf_string != None:
                f.write(">" + append_to_scaf_string + str(key))
                f.write('\n')
                f.write(str(val.seq))
                f.write('\n')
            elif append_to_scaf_string == None:
                f.write(">" + str(key))
                f.write('\n')
                f.write(str(val.seq))
                f.write('\n')
            else:
                pass

def flaten_fasta(fasta_dict, outputname):
    """Often Times when you're playing around with other peoples Fasta, these
    things are wrapped and it makes my life a pain in the ass. So, with this I
    just wanted to write something that will simply unwrap the sequence and
    print it out.

    :fasta: TODO
    :returns: TODO

    """
    try:
        os.remove(outputname)
    except OSError:
        pass
    
    with open(outputname, 'a+') as f:
        for key, val in fasta_dict.items():
            add_carrot = '>' + str(key)
            f.write(add_carrot)
            f.write('\n')
            f.write(str(val.seq))
            f.write('\n')

def get_parser():
    parser = argparse.ArgumentParser(description='Software to read in fasta \
            file and do basic functionality that is often required  ')
    parser.add_argument('-f','--fasta', help='fasta file to read', \
            required=True, dest='f')
    parser.add_argument('-s','--scaffold', help='scaffold to split',\
            required=False,dest='s') 
    parser.add_argument('-l','--locsplit', help='location in ccaffold to split',\
            required=False,dest='l') 

    parser.add_argument('-slist','--scaffoldlist', help='scaffold list to be \
            manipulated Ns about ',\
            required=False,dest='slist') 

    parser.add_argument('-fl1','--filter1', help='filter1. Sequences that do\
            not meet this cut off will be filtered out. This is the LOW end',\
            required=False,dest='fl1') 

    parser.add_argument('-fl2','--filter2', help='filter2. Sequences that are\
            longer than this value will be sequenced out',\
            required=False,dest='fl2') 

    parser.add_argument('-ra','--renameappend', help='append string to\
            scaffolds to rename',
            required=False,dest='ra') 
    parser.add_argument('-uw','--unwrap', help='unwrap fasta file',
            required=False,dest='uw') 
    parser.add_argument('-p','--pull', help='pull out specific scaffolds, write \
            to file',required=False,dest='p') 

    parser.add_argument('-st','--stats', help='Reports basic stats mean, \
            meadien, largest, smallest scaf ',required=False, \
            action='store_true', dest='st') 
    
    parser.add_argument('-plt','--pltn', help='Plots Ns on each scaffold,' \
            ,required=False, action='store_true', dest='plt') 

    parser.add_argument('-o','--output', help='output file to write to', \
            required=False, dest='o')
    args = vars(parser.parse_args())    
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    StartTime = datetime.now()
    
    #Read Fasta File
    file_dict = read_in_fasta(args.f)
    
    ##At some point these need to be fixed
    if args.s != None and args.l == None:
        #exit if split given without where
        print("If splitting file need location so split scaffold")
        sys.exit(-2)
    elif args.s == None and args.l != None:
        print("If splitting file need scaffold to split")
        sys.exit(-2)

    #Break scaffold with many different options
    elif args.ra != None and args.s != None and args.l != None:
        broken_scaffold = split_scaffold(file_dict, args.s, args.l)
        brokenscafwriter(args.o,file_dict,args.s,broken_scaffold,args.ra)
    elif args.ra != None and args.s == None and args.l == None:
        rename_scaf_writer(args.o,file_dict,args.ra)
    elif args.ra == None and args.s != None and args.l != None:
        broken_scaffold = split_scaffold(file_dict, args.s, args.l)
        brokenscafwriter(args.o,file_dict,args.s,broken_scaffold,None)

    #If we just want to unwrap fasta
    elif args.uw != None:
        flaten_fasta(file_dict,args.o)
    
    #Filter scaffolds of certain length
    elif args.fl1 != None and args.fl2 != None:
        final_list = filter_fasta_dict(file_dict,args.fl1, args.fl2)
        flaten_fasta(final_list, args.o)
    
    #Stat culculatsion 
    elif args.st != False:
        report_fasta_stats(file_dict)
    
    elif args.plt != None:
        from plt_scaf import *
        chrom_len_dict = calculate_chrom_len(file_dict)
        ns_window_dict = calc_ns_in_window(file_dict, 100000)

        if args.slist != None:
            scaffold_list = args.slist.split(',')
            hist_plot_ns(chrom_len_dict,ns_window_dict,100000,scaffold_list)
        else:
            hist_plot_ns(chrom_len_dict,ns_window_dict, 100000, None)



   
endtime = datetime.now()
finaltime = endtime - StartTime 

print ("Total Time %s" % (finaltime))

