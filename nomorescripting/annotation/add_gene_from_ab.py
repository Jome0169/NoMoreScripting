#!/usr/local/bin/python3
import getopt
import sys
import os
import re
import copy
from sys import argv
from datetime import datetime
import argparse


def read_abinit_names(file_abinit):
    """reads in abinit names into a list

    :file_abinit: TODO
    :returns: TODO

    """
    abinit_names = []
    with open(file_abinit, 'r') as f:
        for line in f:
            clean = line.strip()
            abinit_names.append(clean)
    return abinit_names

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


def get_true_name(ID_string, num):
    """TODO: Docstring for get_true_name.

    :arg1: TODO
    :returns: TODO

    """
    if num == 1:
        split_ID_string = ID_string.split(';')
        mrna_name_reformat = copy.deepcopy(split_ID_string[1]).replace('Name=','')
        return mrna_name_reformat
    elif num == 2:
        split_ID_string = ID_string.split(';')
        mrna_name_reformat =copy.deepcopy(split_ID_string[2]).replace('Target=','')
        return mrna_name_reformat



def reformat_gff_id(ID_string, part):
    """TODO: Docstring for reformat_gff_id.

    :ID_string: TODO
    :part: TODO
    :new_name: TODO
    :returns: TODO

    """
    
    if part == 'gene':

        split_ID_string = ID_string.split(';')
        name_reformat = copy.deepcopy(split_ID_string[1]).replace('-mRNA-1','').replace('abinit','processed').replace('Name=','')
        final_string_reformat = 'ID=' + name_reformat +  ';' +'Name=' + name_reformat
        return final_string_reformat
    
    elif part == 'mRNA':

        split_ID_string = ID_string.split(';')
        parent_name_reformat = copy.deepcopy(split_ID_string[1]).replace('-mRNA-1','').replace('abinit','processed').replace('Name=','') 
        mrna_name_reformat = copy.deepcopy(split_ID_string[1]).replace('abinit', 'processed').replace('Name=','')

        final_string_reformat = 'ID=' + mrna_name_reformat + ';' + 'Parent=' + \
        parent_name_reformat + ';' +  'Name=' + mrna_name_reformat +';' + ';'.join(split_ID_string[2:])
        return final_string_reformat

    elif part == 'CDS':

        split_ID_string = ID_string.split(';')
        mrna_name_reformat = copy.deepcopy(split_ID_string[2]).replace('abinit','processed').replace('Target=','')
        final_string_reformat = 'ID=' + mrna_name_reformat + ':' + "cds" + ';' + 'Parent=' + mrna_name_reformat
        return final_string_reformat

def fix_gff_names(abinit_list, gff_list, output_file):
    """takes in both the gff and abinit list and replaces poorly annoated
    abinit gene names with the markings of a gene in gff3 format

    :abinit_list: TODO
    :gff_list: ['Chr6', 'maker', 'exon', '131200', '131899', '.', '+', '.',
    'ID=maker-Chr6-snap-gene-0.137-mRNA-1:exon:31;Parent=maker-Chr6-snap-gene-0.137-mRNA-1']
    :returns: fixed_gff_list 

    """

    #fixed_gff_list = []


    remove_files(output_file)

    with open(output_file, 'a+') as f:

        for item in gff_list:
            if len(item) > 2:
                
                if item[2] == "match":
                    #Check name function looks and isolates the mrna anme nad
                    #checks it against list 
                    check_name = get_true_name(item[8], 1)
                    if check_name in abinit_list:
    
    
                        new_gene_format = copy.deepcopy(item)
                        new_gene_format[1] = 'maker'
                        new_gene_format[2] = 'gene'
                        gene_fixed_line  = reformat_gff_id(new_gene_format[8],'gene')
                        new_gene_format[8] = gene_fixed_line
    
                        new_mrna_format = copy.deepcopy(item)
                        new_mrna_format[1] = 'maker'
                        new_mrna_format[2] = 'mRNA'
                        mrna_fixed_line = reformat_gff_id(new_mrna_format[8],'mRNA')
                        new_mrna_format[8] = mrna_fixed_line
    
                        #bad_good_good_list = [item, new_gene_format,new_mrna_format]
                        #fixed_gff_list.append(bad_good_good_list)
    
                        f.write(('\t'.join(new_gene_format)))
                        f.write('\n')
                        f.write(('\t'.join(new_mrna_format)))
                        f.write('\n')

                elif item[2] == "match_part":
                    check_name = get_true_name(item[8], 2)
                    
                    if check_name in abinit_list:
    
                        new_cds_format = copy.deepcopy(item)
                        new_cds_format[1] = 'maker'
                        new_cds_format[2] = 'CDS'
                        CDS_fixed_line = reformat_gff_id(new_cds_format[8],'CDS')
                        new_cds_format[8] =CDS_fixed_line
    
                        #bad_good_good_list = [item, new_cds_format]
                        #fixed_gff_list.append(bad_good_good_list)
                         
                        f.write(('\t'.join(new_cds_format[0:9])))
                        f.write('\n')

                else:
                    quick_fix = ' '.join(item[8:])
                    item[8] = quick_fix
                    f.write(('\t'.join(item[0:9])))
                    f.write('\n')
            else:
                if len(item) == 2:
                    f.write(' '.join(item))
                    f.write('\n')
                else:
                    f.write(item[0])
                    f.write('\n')




def remove_files(file_name):
    """removes file name if it exits

    :file_name): TODO
    :returns: TODO

    """

    try:
        os.remove(file_name)
    except OSError:
        pass




def get_parser():
    parser = argparse.ArgumentParser(description='The purpose of this script is \
            add in abinitio gene predictions back into the mix in AllMaps.')
    parser.add_argument('-g','--gff', help='gff file to add abinitio genes back into', required=True, dest='g')
    parser.add_argument('-n','--names', help='names of abinit genes', required=True, dest='n')

    parser.add_argument('-o','--output', help='File to Write output to. Will \
            Remove file if finds a conflicting name', required=True, dest='o')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    StartTime = datetime.now()
    args = get_parser().parse_args()

    gene_names = read_abinit_names(args.n)
    gff_names = read_in_gff(args.g)

    fix_gff_names(gene_names,gff_names, args.o)



#Speed Things
EndTime = datetime.now()
FinalTime = EndTime - StartTime

print ("Total Time %s" % (FinalTime))



