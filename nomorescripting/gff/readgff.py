from datetime import datetime
from sys import argv
import copy
import argparse
import sys
import os
import itertools



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

def retrieve_max_scaf_len(gff_record):
    """returns a nested list of scaffodls with the largest number found in
    scaffold. Will work as a false representation when splitting gff scaffold

    :gff_record: TODO
    :returns: TODO

    """
    scaffold_list = []
    for item in gff_record:
        if item[0] not in scaffold_list:
            scaffold_list.append(item[0])
    
    nested_scaf_list = [] 
    for scaf in scaffold_list:
        nested_scaf_list.append([scaf])
    
    for hit in gff_record:
        for scaf_name in nested_scaf_list:
            if hit[0] in scaf_name:
                scaf_name.append(int(hit[3]))
                scaf_name.append(int(hit[4]))

    final_scaf_list = []
    for finaliter in nested_scaf_list:
        largesthit = max(finaliter[0:])
        final_scaf_list.append(finaliter[0])
        final_scaf_list.append(largesthit)
    return final_scaf_list 

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


def write_output(gff_to_write, output_file):
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
        write_output(fixed_gff_file,args.o)

    elif args.s != None and args.l != None and args.ra != None:
        fixed_gff_file = fix_split_scaffold(gff_file,args.s, args.l)
        rename_scaffold(fixed_gff_file,args.ra)
        write_output(fixed_gff_file,args.o)









