
"""
Takes in an AHRD output file, and a gff file with corresponding names and adds
the function from the AHRD file into the Gff3 file.

"""
from sys import argv 
from datetime import datetime
from read_gff import read_in_gff
from sys import argv
import copy
import argparse
import sys
import os


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
                annotation_key[clean_line[0]] = clean_line[-1]

    return annotation_key

def read_mapping_file(map_file):
    """reads in mapping file

    :map_file: TODO
    :returns: TODO

    """
    new_name_old_name = {}
    with open(map_file, 'r') as z:
        for line in z:
            clean = line.strip().split()

            new_name_old_name[clean[1]] = clean[0]

    return new_name_old_name


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


def edit_gff_list(gff_list, ahrd_dict):
    """TODO: Docstring for edit_gff_list.

    :gff_list: TODO
    :map_dict: CsGy6G000010 snap_masked-Chr6-processed-gene-0.58
               CsGy6G000020.1 snap_masked-Chr6-processed-gene-0.73-mRNA-1
    :fasta_dict: sp|Q197D5|025R_IIV3 :  Uncharacterized protein 026R OS=Invertebrate iridescent virus 6 OX=176652 GN=IIV6-026R PE=4 SV=1
    :blast_dict: CsGy6G000110.1 : sp|Q8RWB8|UPL6_ARATH
    :returns: TODO

    """

    for gff in gff_list:
        #Check for ## signs denotating new gene
        if len(gff) > 2:
            if gff[2] == 'gene':

                reformat_ID_string_base = gff[8]
                key_name = reformat_ID_string_base.split(';')[0].replace('ID=', '') + '.1'

                clean_out_AED = reformat_ID_string_base.split(';')
                clean_ID_Name=(';'.join(clean_out_AED[0:2]))
                
                if key_name in ahrd_dict:
                    annoation_hit = ahrd_dict[key_name].replace('=', '').replace(';','')

                    create_funct_string = ";Note=%s" % annoation_hit 
                    final_func_string = clean_ID_Name + \
                    create_funct_string

                    gff[8] = final_func_string
                    print('\t'.join(gff))
                else:
                    pass
                    print('\t'.join(gff))

            elif gff[2] == 'mRNA':

                reformat_ID_string_base = gff[8]
                key_name = reformat_ID_string_base.split(';')[0].replace('ID=', '')

                clean_out_AED = reformat_ID_string_base.split(';')
                #Conserve ID parent Name
                clean_ID_Name=(';'.join(clean_out_AED[0:3]))
 

                if key_name in ahrd_dict:
                    annoation_hit = ahrd_dict[key_name].replace('=', '').replace(';','')
                    
                    create_funct_string = ";Note=%s" % annoation_hit 

                    final_func_string = clean_ID_Name + \
                    create_funct_string
                    gff[8] = final_func_string
                    print('\t'.join(gff))

                else:
                    print('\t'.join(gff))


            elif gff[0].startswith('#'):
                print('\t'.join(gff))


                    
            elif gff[2] != 'gene' or gff[2] != 'mRNA' and '#' not in gff[0]:
                safe_copy = gff[8]
                quick_reformat = ''.join(safe_copy)
                gff[8] = quick_reformat
                print('\t'.join(gff))


            else:
                #Sometimes the non gene regions have weird formatting on the end
                reformat_list = gff[8:]
                ' '.join(reformat_list)

                del gff [8:]

                gff[8] = reformat_list
                print('\t'.join(gff))

        else:

            print(' '.join(gff))


def get_parser():
    parser = argparse.ArgumentParser(description='Adds functional annotation to \
            gff3 file. Requires AHRD output in csv format, and a gff3 file with \
            corresponding names')
    parser.add_argument('-g','--gff', help='gff file to read', \
            required=True, dest='g')
    parser.add_argument('-ahrd','--ahrd', help='prefix to ',\
            required=True,dest='ahrd') 
    args = vars(parser.parse_args())    
    return parser




if __name__ == "__main__":
    args = get_parser().parse_args()
    annotation_key = read_ahrd_csv
    read_gff_file = read_in_gff(args.g)
    edit_gff_list(gff_lists,annotation_key)






























