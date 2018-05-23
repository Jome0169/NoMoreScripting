

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



def hist_gene_found(hist_list, output_file):
    """Takes in melted dictionary frmo calculate_gene_scores, and then create a
    histogram with the percentge of genes found. BEcause we know we've revoeced
    many genes, however, what quality are they?

    :melted_dict: [gene_total_len,percentage_found,final_percent_ID]
    :output_file: TODO
    :returns: TODO

    """
    logging.info("Plotting List")
    bin_number = 100
    


    fig, ax = plt.subplots()    

    ax.set_title('Distribtuion of total length found in %s genes' \
            %len(hist_list))

    n, bins, patches = ax.hist(hist_list, bin_number, density=1)
    plt.show()




def get_parser():
    parser = argparse.ArgumentParser(description="Script that takes in a PAF \
            output format from Minimap2 and returns a double nested dictionary \
            where the gene name is the key, the scaffold name is the second \
            key, and their corresponding hits are the values in a list")
    
    #Arguments to Parse
    parser.add_argument('-mm', '--minimap', help="Minimap 2 output PAF", 
            required=True, dest='mm')

    return parser



if __name__ == "__main__":
    
    args = get_parser().parse_args()
    StartTime = datetime.now()
    
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    
    logging.info('About to read PAF file')
    raw_PAF_file = parse_PAF_file(args.mm) 
    nested_paf_file = calculate_gene_scores(raw_PAF_file)

    gene_len_found = []
    percent_ID = []

    for key, val in nested_paf_file.items():
        for scaffold, hit in val.items():
            gene_len_found.append(float(hit[1]))
            percent_ID.append(float(hit[2]))

    
    #Plot the by
    hist_gene_found(gene_len_found,None)
    hist_gene_found(percent_ID,None)


