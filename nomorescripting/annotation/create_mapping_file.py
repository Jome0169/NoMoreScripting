# -*- coding: utf-8 -*-
"""
    annotation.create_mapping_file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This software creates a maping file from a gff file that can then be used
    for efficient functional annotations. The mapping file will follow
    standared TAIR11 practices. Before running this software ensure that the
    gff fed into this has been SORTED and ORDEERED, as gene names are entiryly
    dependent on locaiton

    :copyright: (c) 2018 by YOUR_NAME.
    :license: LICENSE_NAME, see LICENSE for more details.



   NOTE  Wed Jun 20 16:22:02 EDT 2018: THis thing needs to be edited so it cna
   work with NON chromsome level assembies, as well as increase its general
   flexability
"""

from datetime import datetime
from readgff import read_in_gff
from sys import argv
import copy
import argparse
import sys
import os

def extract_gene_mrna_names(gff_list, prefix_name, output_file):
    """reads in a list of gff3 lines, takes list and parses out either mrNAs or
    the Gene names. After this creates a new name, puts them in seperate list

    :gff_list: TODO
    :returns: TODO

    """
    BaseName = str(prefix_name)

    mRNACounter = 10
    geneCounter = 10

    total_name_count = 6
    ChromNumber = None

    for item in gff_list:
        
        if len(item) == 2:
            pass
        else:
            if item[2] == 'mRNA':
                
                take_chrom_num = item[0].replace('Chr', '')
                if ChromNumber == None:
                    ChromNumber = take_chrom_num
                elif ChromNumber != None and take_chrom_num != ChromNumber:
                    ChromNumber = take_chrom_num
                    mRNACounter = 10
                    geneCounter = 10
                elif ChromNumber != None and take_chrom_num == ChromNumber:
                    pass

                take_chrom_num = item[0].replace('Chr', '')

                clearsplit = item[8].replace(';',':').split(':')
                real_name = clearsplit[0].replace('ID=', "")


                new_format = BaseName + str(take_chrom_num) + 'G' + str(mRNACounter).zfill(int(6)) + '.1'

                old_name_raw = real_name.split('=')


                if output_file == None: 
                    print(old_name_raw[0],'\t',new_format)
                elif output_file != None:
                    with open(output_file, 'a+') as f:
                        create_final_string = str(old_name_raw[0] + '\t' + new_format)
                        f.write(create_final_string)
                        f.write('\n')

                mRNACounter += 10



            elif item[2] == 'gene':

                take_chrom_num = item[0].replace('Chr', '')
                if ChromNumber == None:
                    ChromNumber = take_chrom_num
                elif ChromNumber != None and take_chrom_num != ChromNumber:
                    ChromNumber = take_chrom_num
                    mRNACounter = 10
                    geneCounter = 10
                elif ChromNumber != None and take_chrom_num == ChromNumber:
                    pass

                take_chrom_num = item[0].replace('Chr', '')

                clearsplit = item[8].replace(';',':').split(':')
                real_name = clearsplit[1].replace('ID=', "")

                new_format = BaseName + str(take_chrom_num) + 'G' + str(mRNACounter).zfill(int(6))
                
                old_name_raw = real_name.split('=')

                if output_file == None: 
                    print(old_name_raw[1],'\t',new_format)
                elif output_file != None:
                    with open(output_file, 'a+') as f:
                        create_final_string = str(old_name_raw[1] + '\t' + new_format)
                        f.write(create_final_string)
                        f.write('\n')


                geneCounter += 10


def get_parser():
    parser = argparse.ArgumentParser(description='Reads in gff gene file and \
            creates mapping file for rename. This is to be used downstream with\
            other script that were designed to rename all files. Note: This \
            program needs a FINALIZED genome with chromosome numbers and not \
            scaffods. That is the only way a correct mapping file will be\
            generated.')
    parser.add_argument('-g','--gff', help='gff file to read', \
            required=True, dest='g')
    parser.add_argument('-pre','--prefix', help='prefix to ',\
            required=True,dest='pre') 
    parser.add_argument('-o','--output', help='output file to write to', \
            required=False, dest='o')
    args = vars(parser.parse_args())    
    return parser



if __name__ == "__main__":
    args = get_parser().parse_args()
    read_gff_file = read_in_gff(args.g)
    if args.o == None:
        extract_gene_mrna_names(read_gff_file, args.pre, None)
    else:
        extract_gene_mrna_names(read_gff_file, args.pre, args.o)

