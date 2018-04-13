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


def read_extra_annotate_file(prot_blast_inf):
    """TODO: Docstring for read_extra_annotate_file.

    :arg1: TODO
    :returns: TODO

    """
    
    annotations_to_add = {}
    with open(prot_blast_inf, 'r') as f:
        for line in f:
            clean_line = line.strip()
            if clean_line.startswith("CsGy"):
                cucumber_name_digest = clean_line.split(' ')
                key_name = cucumber_name_digest[0]
                annotations_to_add[key_name] = []
            elif clean_line.startswith('Gy14'):
                annotations_to_add[key_name].append(clean_line.split())
            else:
                pass
    return annotations_to_add 


def read_arabidopsis_fasta(arab_file):
    """reads arabidopsis headers into 

    :arab_file: TODO
    :returns: TODO

    """
    arab_func_Dict = {}
    with open(arab_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                clean = line.strip().replace('>','').split('|')
                arab_func_Dict[clean[0]] = clean[1]
    return arab_func_Dict            


def get_arab_info(arab_dict,hitname):
    """TODO: Docstring for get_arab_info.

    :arab_dict: TODO
    :hitname): TODO
    :returns: TODO

    """
    functional_string = ""
    for gene_id, function in arab_dict.items():
        if hitname in gene_id:
            functional_string = function
            break
        else: 
            functional_string = "Uncharacterized protein"
    return functional_string


def filter_prot_file(extra_prot_info, arabidop_dict):
    """TODO: Docstring for coombine_info.

    :ahrd_dict: TODO
    :extra_prot_info: TODO
    :arabidop_dict: TODO
    :returns: TODO

    """

    final_prot_info_dict = {}

    for genename, blast_results in extra_prot_info.items():
        if len(blast_results) == 0:
            final_prot_info_dict[genename] = 'Unknown Protein'
        elif len(blast_results) != 0:
            for item in blast_results:
                if item[1].startswith('AT'):
                    possible_funct = get_arab_info(arabidop_dict,item[1])
                    final_prot_info_dict[genename] = possible_funct
                    break
                else:
                    final_prot_info_dict[genename] = "Uncharacterized protein"
    return final_prot_info_dict


def fix_arabidop_AHRD(ahrd_dict, arab_dict):
    """TODO: Docstring for fix_arabidop_AHRD.

    :ahrd_dict: TODO
    :returns: TODO

    """
    for key, val in ahrd_dict.items():
        if len(val) > 0:
            if val[0].startswith('AT'):
                arab_hit = val[0]
                real_funct = get_arab_info(arab_dict, arab_hit)
                val[-1] = real_funct 
    return ahrd_dict




def edit_ahrd_file_final(ahrd_dict,final_protein_info):
    """TODO: Docstring for edit_ahrd_file_final.

    :ahrd_dict: TODO
    :final_protein_info: TODO
    :returns: TODO

    """
    for key, val in ahrd_dict.items():
        if key in final_protein_info:
            new_info = final_protein_info[key]
            if new_info != val[-1]:
                val[-1] = new_info
    return ahrd_dict

def write_output(final_edited_ahrd,output_name):
    """writes output file. Will remove output if found.

    :output_name): TODO
    :returns: TODO
    """
    remove_files(output_name)
    
    with open(output_name, 'a') as f:
        for key, val in final_edited_ahrd.items():
            val.insert(0,key)
            f.write('\t'.join(val))
            f.write('\n')


def remove_files(file_name):
    """Removes function quickly and easily.

    """
    try:
        os.remove(file_name)
    except OSError:
        pass



def get_parser():
    parser = argparse.ArgumentParser(description='The purpose of this script iA \
            to add additional inforamtion into the ahrd_output file. There are\
            issues with naming proteins "Uncatagorized" vs Unknown. This script\
            will add back in Uncategorized using the input from another script\
            Filter.script.sh.')
    parser.add_argument('-a','--ahrd', help='AHRD output file', required=True, dest='a')
    parser.add_argument('-p','--plant', help='arabiposis prot fasta seq', required=true, dest='p')
    parser.add_argument('-ip','--interpro', help='interproscan file to add to output', required=true, dest='ip')
    parser.add_argument('-go','--goterms', help='go term output file to add to AHRD output file', required=true, dest='go')

    parser.add_argument('-i','--info', help='Takes in the check.prot.info\
            output from bash script filter.script.sh', required=True, dest='i')


    parser.add_argument('-o','--output', help='File to Write output to. Will \
            Remove file if finds a conflicting name', required=True, dest='o')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    StartTime = datetime.now()
    args = get_parser().parse_args()

    ahrd_dict = read_ahrd_csv(args.a)
    plant_info = read_arabidopsis_fasta(args.p)
    ammend_proteins = read_extra_annotate_file(args.i)
    fix_arabidop_AHRD(ahrd_dict,plant_info)

    final_prot_info_dict = filter_prot_file(ammend_proteins,plant_info)
    edit_ahrd_file_final(ahrd_dict,final_prot_info_dict)

    write_output(ahrd_dict, args.o)

   



#Speed Things
EndTime = datetime.now()
FinalTime = EndTime - StartTime

print ("Total Time %s" % (FinalTime))


