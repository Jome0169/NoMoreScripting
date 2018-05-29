# -*- coding: utf-8 -*-
"""
    alignment.PAF_analysis
    ~~~~~~~~~~~~~~~~~~~~~~

    This module analyzes the PAF output format provided from minimap2 output.
    I'm currently using this to analyze what percentage of genes were captured
    in the de-novo assembies.

    

    :copyright: (c) 2018 by YOUR_NAME.
    :license: LICENSE_NAME, see LICENSE for more details.
"""





from datetime import datetime
import logging 
import copy
import argparse
import sys
import os
import itertools
import statistics as s
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
import logging
from PAF_reader import *



def hist_from_list(hist_list,title,output):
    """ Analysis of alignments is a huge deal. This small function just
    createas a basic histogram in order to allow me to visualzie the
    distribution of certain values more easily.
    :output_file: Creates a histogram from list of values

    """


    logging.info("Plotting List")
    bin_number = 100

    fig, ax = plt.subplots()    

    ax.set_title('%s  n=%s' \
            % (title,len(hist_list)))

    n, bins, patches = ax.hist(hist_list, bin_number, density=1)


    plt.show()

    if output == None:
        logging.info("No output provided, showing plots")
        plt.show()
    elif output != None:
        final_string = output + '.pdf'
        fig.savefig(final_string)



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
    #hist_from_list(len_dist,'Number of times each gene was found', "dist.times")
    return(genes_sinlge,genes_more_than_once)



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
                final_list = [gene_name, scaffold, group, \
                str(len(hit_info)),str(info[0]), str(info[1]), str(info[2])]
                f.write(','.join(final_list))
                f.write('\n')


def get_parser():
    parser = argparse.ArgumentParser(description="Script that takes in a PAF \
            output format from Minimap2 and returns a double nested dictionary \
            where the gene name is the key, the scaffold name is the second \
            key, and their corresponding hits are the values in a list")
    
    #Arguments to Parse
    parser.add_argument('-mm', '--minimap', help="Minimap 2 output PAF", 
            required=True, dest='mm')

    parser.add_argument('-o', "--ouptut", help="Output basename to write file\
            to", required=False, dest='o')

    return parser



if __name__ == "__main__":
    
    args = get_parser().parse_args()
    StartTime = datetime.now()
    
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    
    logging.info('About to read PAF file')
    raw_PAF_file = parse_PAF_file(args.mm) 

    #mapping_file_dist = []
    #for hit, locaion in raw_PAF_file.items():
    #    for scaffold, hitinfo in locaion[0].items():
    #        for thing in hitinfo:
    #            mapping_file_dist.append(int(thing[11]))
    #hist_from_list(mapping_file_dist, "Dist of mapping Quals", None)

    #Filter mapping scores, melt dict
    nested_paf_file = calculate_gene_scores(raw_PAF_file)

    single_gene_dict,mutli_gene_dict = seperate_into_classes(nested_paf_file)
    
    write_output_table(mutli_gene_dict,"multi", args.o)    
    write_output_table(single_gene_dict,"single", args.o)    

   
