from datetime import datetime
import copy
import argparse
import sys
import os
import itertools
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py
import logging
#Hacky as heck
from PAF_analysis import hist_from_list

def only_mRNA(gff3_list):
    """Takes in nested gff3 list from readgff, and only returns list with
    'mrna' in the third position. This being due to the fact that gmap stores
    all of it's quality information in the ID column within the gff3 file.
    These filtered mRNa will later be converted into a dict object

    :gff3_list: TODO
    :returns: TODO

    """
    mRNA_items = []

    for item in gff3_list:
        if item[2] == 'mRNA':
            mRNA_items.append(item)

    return mRNA_items



###Think about this function more.
def clobber_mRNA_list(mRNA_list):
    """Takes in the mRNA list, and digest it so the return is a key value
    format where the gene will be the key, and the value will be other gff3
    hits.

    :mRNA_list: TODO
    :returns: TODO

    """


def print_cheker(thing_print):
    """I spend alot of time doing print check on things. This is a way to try
    and automate all the thins I end up having to print check. 

    :thing_print: TODO
    :returns: TODO

    """

    ID_type = type(thing_print)

    if ID_type == "list":
        print(thing_print)




def read_blat_PSL(PSL_file):
    """Reads in PSL file. Reads this in to a list object, and returns this list
    for later processing. 

    :PSL_file: TODO
    :returns: TODO

    """
    
    PSL_list = []
    try:
        open(PSL_file, 'r')
    except ValueError:
        print("Error! No %s PSL to read. Does this file exists?" % PSL_file)

    with open(PSL_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split()
            PSL_list.append(clean_line)
    return PSL_list

def create_gene_dict(PSL_list):
    """Takes in the PSL list, and creates a dictionary of all of the query
    values. This dictionary is to be later used in order to generate a nested
    dictionary where the final format will be

    {query: {target1:[hit1, hit2, hit3], target2:[hit1, hit2, hit3], etc...}}

    :PSL_list: TODO
    :returns: TODO

    """
    query_dict = {}

    for item in PSL_list:
        query_name = item[9]
        query_dict[query_name] = {}

    return query_dict

def nest_target_gene_dict(query_dict, PSL_list):
    """Iterates through the PSL list again, and creates a nested dictionary
    object.

    :query_dict: TODO
    :PSL_list: TODO
    :returns: {query: {target1:[hit1, hit2, hit3], target2:[hit1, hit2, hit3], etc...}}

    """
    #This way of creating a nested dictionary isn't actually garbage. Pretty
    #slick 

    for item in PSL_list:
        query_name = item[9]
        target_name = item[13]
    
        inner_query_dict = query_dict[query_name]
        
        if target_name not in inner_query_dict:
            inner_query_dict[target_name] = [item]
        elif target_name in inner_query_dict:
            inner_query_dict[target_name].append(item)

    return query_dict 

def clobber_inner_list(neseted_query_dict):
    """Takes in nested PSL dictionary, and clobbers the most inner list into a
    string that can then actualyy be used to graph the overall amount of CDS
    retrieved. This has two advantages. Namely we don't have to parse this ugly
    ass innter list later, and it will make calculating overall CDS quality
    much easier.

    input_inner_list: ['2086', '5', '0', '0', '1', '6', '1', '264', '-', 'gala03g06192', '2118', '0', '2097', 'NODE_29262_length_3826_cov_23.417445', '3826', '1236', '3591', '3', '33,1760,298,', '21,54,1820,', '1236,1533,3293,']]
    :neseted_: {query: {target1:[hit1, hit2, hit3], target2:[hit1, hit2, hit3], etc...}}

    :returns: Something nice

    """
    #copying and editing dict gets werid, so make a deepcopy
    copy_dict = copy.deepcopy(neseted_query_dict)  
    
    counter = 0
    other = 0
    for gene, nested_dict in copy_dict.items():
        for scaffold, hit_info in nested_dict.items():
            reformulat_nested_list = []
            for item in hit_info:

                #Assign vals from PAL
                number_mismatches = int(item[1])
                length_of_query = int(item[10])
                start_query = int(item[11])
                stop_query = int(item[12])
            
                #Calc stats
                alignment_length = (stop_query - start_query)
                percent_ID = (alignment_length - number_mismatches)/length_of_query

                reformulate_list = [length_of_query,alignment_length, percent_ID, ]
                reformulat_nested_list.append(reformulate_list)
            
            nested_dict[scaffold] = reformulat_nested_list
    return copy_dict

def seperate_into_classes(melted_dict):
    """Takes in melted dict, and created two seperate list based off how many
    scaffolds each gene hits to. I am curious to see how many genes hit
    multiple times, as well as how many genes hit once. Should be rather
    interesting.

    :melted_dict: TODO
    :returns: TODO

    """
    len_dist = []
    genes_sinlge = {}
    genes_more_than_once = {}
    
    for key,val in melted_dict.items():
        take_gene_count = len(val)
        len_dist.append(take_gene_count)

        if take_gene_count == 1:
            genes_sinlge[key] = val
        elif take_gene_count > 1:
            genes_more_than_once[key] = val
    
    #plot genes    
    hist_from_list(len_dist,'Number of times each gene was found', "dist.times")
    return(genes_sinlge,genes_more_than_once)


def plot_percent_id(gene_dict):
    """Takes in gene dictiionary, creates a list of all percent ID's and plots
    them. Easy peasy.
    
    [length_of_query,alignment_length, percent_ID, ]
    :gene_dict: TODO
    :returns: TODO

    """
        
    gene_PID_list = []
    print(len(gene_dict))
    for gene, hits in gene_dict.items():
        for scaffold, info in hits.items():
            for item in info:
                gene_PID_list.append(item[2])

    hist_from_list(gene_PID_list, "Percent ID distribution", None)


def write_output_table(gene_dict, group, output_name):
    """Writes the nested hit information into a CSV table output for ease of
    analysis and computation. Scripting makes data analysis far more
    challenging than it needs to be, so this ensures ease of analysis and
    compariong between two different alginers. GMAP, and minimap

    :gene_dict: Nested dict with hits
    :group: single or multi
    :output_name: output base file name
    :returns: TODO

    """
    with open(output_name, 'a') as f:
        for gene_name, hit_info in gene_dict.items():

            for scaffold, info in hit_info.items():
                for i in info:

                    final_list = [gene_name, scaffold, group, str(len(hit_info)),
                            str(i[0]), str(i[1]), str(i[2])]
                    f.write(','.join(final_list))
                    f.write('\n')
                    print(final_list)




def get_parser():
    parser = argparse.ArgumentParser(description="Softare that reads in gmap\
            file, and reports back some basic statistics about it. This files\
            is meant to be used with various de novo assemlbies that have CDS\
            sequenced back to them using GMAP. This way we can access the\
            overall quality of a genomic assembly. File requires both gff3 and\
            sam output in 'gene' format as GMAP will not print mapping files\
            otherwise")

    parser.add_argument("-p",'--psl', help="Gff3 file from GMAP \
to be input into file.", required=True, dest='p')

    parser.add_argument("-o",'--output', help="output file to write gene \
            information to.", required=True, dest='o')

    args = vars(parser.parse_args())
    return parser


if __name__ == '__main__':
    args = get_parser().parse_args()

    PSL_list = read_blat_PSL(args.p)
    query_dict =  create_gene_dict(PSL_list)
    
    #Edits given list
    nest_target_gene_dict(query_dict, PSL_list)
    
    #Digest inner list 
    for key,val in query_dict.items():
        if len(val) >= 2:
            print(key)
            for thing, other in val.items():
                print(thing)
                for i in other:
                    print('\t'.join(i))
            input()



    melted_dict = clobber_inner_list(query_dict)

    single_gene_dict,mutli_gene_dict = seperate_into_classes(melted_dict)

    
    



    #write_output_table(single_gene_dict, "single", args.o)
    #write_output_table(mutli_gene_dict, "multi", args.o)

    #gotta fix these plots, and clean them up
    #plot_percent_id(single_gene_dict)
    #plot_percent_id(mutli_gene_dict)
