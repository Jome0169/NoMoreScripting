from datetime import datetime
from sys import argv
import copy
import argparse
import sys
import os
import itertools
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py


def read_in_gff(arg1):
    """TODO: Docstring for read_in_gff.

    :arg1: TODO
    :returns: TODO

    """
    all_gff_records = []
    
    with open(arg1, 'r') as f:
        for line in f:
            if '#' in line:
                pass
            else:
                cleanline = line.strip().split()
                all_gff_records.append(cleanline)
    return all_gff_records

def read_in_gff2(arg1):
    """TODO: Docstring for read_in_gff.

    :arg1: TODO
    :returns: TODO

    """
    all_gff_records = []
    
    with open(arg1, 'r') as f:
        parse_smart = []
        for line in f:
            if "##gff-version" in line:
                pass
            elif "##gff-version" not in line and '#' not in line:
                cleanline = line.strip().split()
                parse_smart.append(cleanline)

            elif "##gff-version" not in line and '#' in line:
                all_gff_records.append(parse_smart)
                parse_smart = []
            elif line == '' and parse_smart != None:
                all_gff_records.append(parse_smart)
            elif line == '' and parse_smart == None:
                pass
   
    #Remove Empties
    no_empty_list = [x for x in all_gff_records if x != []]
    return no_empty_list


def rename_scaffold(gff_record, scaf_rename):
    """TODO: Docstring for rename_scaffold.

    :gff_record: TODO
    :scaf_rename: TODO
    :returns: TODO

    """
    for item in gff_record:
        reformat_name = item[0]
        new_scaf_name = str(scaf_rename) + reformat_name
        item[0] = new_scaf_name
    return gff_record



def calculate_gene_len(nested_gff_list):
    """TODO: Docstring for calculate_gene_len.

    :gff_record: TODO
    :returns: TODO

    """
    total_gene_len = None

    for annotation_hit in nested_gff_list:
        if annotation_hit[2] == 'mRNA':
            total_len =  int(annotation_hit[4]) - int(annotation_hit[3])
            total_gene_len = total_len
        else:
            pass
    return total_gene_len

def plot_lens(nested_gff_list):
    """TODO: Docstring for plot_lens.

    :something: TODO
    :returns: TODO

    """
    final_list = [] 
    #Get All Lens 
    for item in nested_gff_list:
        z = calculate_gene_len(item)
        final_list.append(int(z))
    #Plot 
    plt.hist(final_list)
    plt.title("Gene Length Histogram")
    plt.xlabel("Gene Length")
    plt.ylabel("Frequency")
    plt.show()


def filter_gff_by_len(nested_gff,smallest,largest):
    """filters out exon hits that do not meet certain length criterias. This is
    useful in that some programs that output gff3 files have very truncated
    gff3 annoations that do not apper real. Takes in three arguemtns, nested
    gff3, smallest value to filter out, largest value to filter out.

    :nested_gff: TODO
    :smallest: TODO
    :largest: TODO
    :returns: TODO

    """
    passing_len_criteria = []
    failing_len_criteria = []

    for item in nested_gff:
        total_gene_len = calculate_gene_len(item)
        if total_gene_len <= int(largest) and total_gene_len >= int(smallest):
            passing_len_criteria.append(item)
        else:
            failing_len_criteria.append(item)
    
    return passing_len_criteria


def fix_split_scaffold(gff_record,split_scaffold,loc_split):
    """TODO: Docstring for fix_split_scaffold.
    :gff_record: TODO
    :split_scaffold: TODO
    :loc_split: TODO
    :returns a reformatted list where gff lines that hit the scaffold that will
    be split have their hit locations recalculated

    """
    reformatted_gff = []
    
    for gff in gff_record:
        if gff[0] != split_scaffold:
            reformatted_gff.append(gff)
        elif gff[0] == split_scaffold:
            #Do Genes hit to certain ones
            if int(gff[3]) < int(loc_split) and int(gff[4]) < int(loc_split):
                rename = gff[0] + '_1'
                gff[0] = rename
                reformatted_gff.append(gff)

            elif int(gff[3]) > int(loc_split) and int(gff[4]) > int(loc_split):
                #If we need to split this scaffold, recalibrate bp position
                rename = gff[0] + '_2'
                gff[0] = rename 
                hit1 = int(gff[3]) - int(loc_split)
                hit2= int(gff[4]) - int(loc_split)
                gff[3] = str(hit1)
                gff[4] = str(hit2)
                reformatted_gff.append(gff)

            elif int(gff[3]) > int(loc_split) and int(gff[4]) < int(loc_split):
                print("Gene in Middle, find alt fix.")
                sys.exit(2)

    return reformatted_gff


def write_output_file(gff_to_write, output_file):
    """writes final list to output file.

    :gff_to_write: TODO
    :output_file: TODO
    :returns: TODO

    """
    try:
        os.remove(output_file)
    except OSError:
        pass

    with open(output_file, 'a+') as f:
        f.write('##gff-version 3')
        f.write('\n')
        for item in gff_to_write:
            f.write('\t'.join(item))
            f.write('\n')


def write_output_file2(gff_to_write, output_file):
    """writes final list to output file.

    :gff_to_write: TODO
    :output_file: TODO
    :returns: TODO

    """
    try:
        os.remove(output_file)
    except OSError:
        pass

    with open(output_file, 'a+') as f:
        f.write('##gff-version 3')
        f.write('\n')
        for nest in gff_to_write:
            for item in nest:
                f.write('\t'.join(item))
                f.write('\n')
            f.write("###")
            f.write('\n')










def get_parser():
    parser = argparse.ArgumentParser(description='Software to read in gff\
            file and do basic functionality that is often required  ')
    parser.add_argument('-g','--gff', help='gff file to read', \
            required=True, dest='g')
    parser.add_argument('-ra','--renameappend', help='Rename base scaffold', \
            required=False, dest='ra')
    parser.add_argument('-s','--scaffold', help='scaffold to update',\
            required=False,dest='s') 
    parser.add_argument('-l','--locsplit', help='location in scaffold that was split',\
            required=False,dest='l') 

    parser.add_argument('-f1','--filter1', help='Remove gff genes smaller than number',\
            required=False,dest='f1') 
    parser.add_argument('-f2','--filter2', help='Remove Gff genes larger than number',\
            required=False,dest='f2') 

    parser.add_argument('-o','--output', help='output file to write to', \
            required=True, dest='o')
    args = vars(parser.parse_args())    
    return parser



if __name__ == "__main__":
    args = get_parser().parse_args()
    gff_file = read_in_gff(args.g)
    
    if args.s == None and args.l != None:
        print("Need Scaffold Name to split")
        sys.exit(2)
    elif args.s != None and args.l == None:
        print("Need Location to split scaffold")
        sys.exit(2)
    #split scaffold no rename
    elif args.s != None and args.l != None and args.ra == None:
        fixed_gff_file = fix_split_scaffold(gff_file,args.s, args.l)
        write_output_file(fixed_gff_file,args.o)

    elif args.s != None and args.l != None and args.ra != None:
        fixed_gff_file = fix_split_scaffold(gff_file,args.s, args.l)
        rename_scaffold(fixed_gff_file,args.ra)
        write_output_file(fixed_gff_file,args.o)
    
    elif args.f1 != None or args.f2 != None:
        tryagain = read_in_gff2(args.g)
        final_list = filter_gff_by_len(tryagain,args.f1,args.f2)
        write_output_file2(final_list,args.o)
        plot_lens(final_list)











