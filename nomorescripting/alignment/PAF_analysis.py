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



def parse_PAF_file(PAF_file):
    """Reads in PAF file output from Minimap2. Minimap2 is used to Align CDS
    sequences to the contiga.fa. This way I can then create a brief analysis of
    how many genes are present and absent, allowing a far better idea of the
    overall quality of the genome.

    :PAF_file: TODO
    :returns: TODO

    """

    def double_nest_dict(PAF_dict):
        """Takes in dictionary with gene name as key and creates a double
        nested dictionary with scaffold being the other internal key. This is
        in the hops of calculating real good gene hits in order Identify how
        well certain de novo assemblies perform with gene capture
        :returns: TODO

        maker-10000558-snap-gene-0.10-mRNA-1 [
        {'NODE_115882_length_1094_cov_3.855199': [['maker-10000558-snap-gene-0.10-mRNA-1', '81', '17', '81', '+', 'NODE_115882_length_1094_cov_3.855199', '1094', '967', '1031', '64', '64', '23', 'tp:A:P', 'cm:i:16', 's1:i:64', 's2:i:50', 'dv:f:0.0079']], 
        'NODE_69614_length_1708_cov_5.135727': [['maker-10000558-snap-gene-0.10-mRNA-1', '81', '17', '78', '+', 'NODE_69614_length_1708_cov_5.135727', '1708', '1581', '1642', '50', '61', '0', 'tp:A:S', 'cm:i:7', 's1:i:50', 'dv:f:0.0592']]}]


        """


        final_PAF_dict = {}
        for key, value in PAF_dict.items():
            #Add Gene Key Back in
            final_PAF_dict[key] = []
            #Create internal Dict 
            create_scaf_dict = {}
            for item in value:
                if item[5] not in create_scaf_dict:
                    create_scaf_dict[item[5]] = [item]
                elif item[5] in create_scaf_dict:
                    create_scaf_dict[item[5]].append(item)
            final_PAF_dict[key].append(create_scaf_dict)

        return final_PAF_dict

    
    try:
        open(PAF_file, 'r')
    except :
        print("PAF file was not given will not report CDS retrivle")
        logging.warning("No CDS given")
        return None

    gene_hit_dict = {}
    with open(PAF_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split('\t')
            if len(clean_line) <= 12:
                #Break if not PAF File
                print("Incorrect ")
                return None
                break

            elif clean_line[0] not in gene_hit_dict:
                gene_hit_dict[clean_line[0]] = [clean_line]
            elif clean_line[0] in gene_hit_dict:
                gene_hit_dict[clean_line[0]].append(clean_line)
            else:
                pass
    
    #Function call to above
    run_this_copy = copy.deepcopy(gene_hit_dict)
    logging.info('Nesting PAF file by gene - scaffold - hits')
    final_nested_dict = double_nest_dict(run_this_copy)
    
    #Return nested dict
    return final_nested_dict



def overlap(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return not (end1 < start2 or end2 < start1)

def calculate_gene_scores(nested_gene_hit_dict):
    """calculate PAF scores for each gene. Calculate each hits length and then add it
    to position 0 in the list. This 

    genemark-10000714-processed-gene-5.52-mRNA-1 [[], 
    ['genemark-10000714-processed-gene-5.52-mRNA-1', '1320', '4', '1317', '+',
    .....'

    :gene_hit_dict: TODO
    :returns: TODO

    maker-10000558-snap-gene-0.10-mRNA-1 {'NODE_115882_length_1094_cov_3.855199': [81, 0.7901234567901234, 1.0], 
    'NODE_69614_length_1708_cov_5.135727': [81, 0.7530864197530864, 0.819672131147541]}
 

    """
    def find_possible_scaf_split(scaffold_lvl_dict):
        """takes in a scaffold level deict and looks for overlap -- or 500 bp
        range of uquery. This is in order to get a better idea of how many CDS
        sequences we're capturing. If we do find this, we create a new
        dictionary addition of ocmbined sets

        :scaffold_lvl_dict: TODO
        :returns: TODO

        """
        merged_split_scaf_dict = {}

        if len(scaffold_lvl_dict) == 2:
            combos = list(itertools.combinations(scaffold_lvl_dict.keys(), 2))
            
            gene1_hits = scaffold_lvl_dict[combos[0][0]][0]
            gene2_hits = scaffold_lvl_dict[combos[0][1]][0]
            
            get_query_location = lambda gene_hit: gene_hit[2:4] 
            
            #Extract Query List aresas
            query_hits_gene_1 = list(map(int,get_query_location(gene1_hits)))
            query_hits_gene_2 = list(map(int,get_query_location(gene2_hits)))
            
            #50 BP range. We're interested in genes that see continououes
            query_hits_gene_1_UPS_DS = [query_hits_gene_1[0]-50,query_hits_gene_1[1] + 50]
            query_hits_gene_2_UPS_DS = [query_hits_gene_2[0]-50,query_hits_gene_2[1] + 50]

            check_overlap = overlap(query_hits_gene_1_UPS_DS[0], query_hits_gene_1_UPS_DS[1], \
                    query_hits_gene_2_UPS_DS[0], query_hits_gene_2_UPS_DS[1])
            
            if check_overlap == True:
                
                new_query_scaffold = combos[0][0] + '____' + combos[0][1]
            
                #Merge list, take smallest and largest
                all_hit_loc = query_hits_gene_1 + query_hits_gene_2 
                start_new_query = min(all_hit_loc)
                end_new_query = max(all_hit_loc)
                
                number_res_matches = int(gene1_hits[9]) + int(gene1_hits[9]) 
                algn_block_len =  ((int(gene1_hits[10]) + int(gene1_hits[10]))/2)

                psuedo_hit = [None] * 12
                
                #Gene Name
                psuedo_hit[0] = gene1_hits[0]
                psuedo_hit[1] = gene1_hits[1]

                #Start Stop sew
                psuedo_hit[2] = str(start_new_query)
                psuedo_hit[3] = str(end_new_query)
                
                #Residue Matches, Algn Block
                psuedo_hit[9] = str(int(number_res_matches))
                psuedo_hit[10] = str(int(algn_block_len))
            
                merged_split_scaf_dict[new_query_scaffold] = [psuedo_hit]

        else:
            #
            for scaf, hit in scaffold_lvl_dict.items():
                merged_split_scaf_dict[scaf] = hit
           
        return merged_split_scaf_dict


    def digest_scaffold_dict(scaffold_lvl_dict):
        """Takes in scaffold dict, and reports back the melted list with the
        alignemtn ID and total len found, etc...

        :arg1: TODO
        :returns: TODO

        maker-10000558-snap-gene-0.10-mRNA-1 {'NODE_115882_length_1094_cov_3.855199': 
        [gene_total_len,percentage_found,final_percent_ID]

        """

        melted_scaffold_hit_dict = {}
        for scaffold_name, hits in scaffold_lvl_dict.items():
            number_hits = len(hits)
            gene_total_len = 0
            total_len_found = 0
            alignment_ID = 0

            for list1 in hits:
                cal_len_hit = int(list1[3]) - int(list1[2])
                alignment_ID_list = int(list1[9]) /int(list1[10])

                total_len_found += cal_len_hit
                gene_total_len = int(list1[1])
                alignment_ID += alignment_ID_list
            
            if gene_total_len == 0 or total_len_found ==0 or alignment_ID ==0:
                pass
            else:
                final_percent_ID = alignment_ID/number_hits 
                percentage_found = total_len_found

                melted_scaffold_hit_dict[scaffold_name] = [gene_total_len,percentage_found,final_percent_ID]

        return melted_scaffold_hit_dict
            
    #Maing function call
    if nested_gene_hit_dict == None:
        return None
    else:
        final_gene_dict = {}
        logging.info("About to melt PAF dict and report back gene_len, \
                percentage_found, and final_percent_ID")
        for gene_hit, value in nested_gene_hit_dict.items():
            #Look for well split genes
            fix_split_scaffs = find_possible_scaf_split(value[0])
            melted_dict = digest_scaffold_dict(fix_split_scaffs)
            final_gene_dict[gene_hit] = melted_dict
    
    logging.info("Finishing dictionary melt")
    return final_gene_dict


def hist_from_list(hist_list,title ):
    """Analysis of alignments is a huge deal. This small function just
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

    #hist_from_list(mapping_file_dist, "Dist of mapping Quals", None)

    #Filter mapping scores, melt dict
    nested_paf_file = calculate_gene_scores(raw_PAF_file)

    #function in this file
    single_gene_dict,mutli_gene_dict = seperate_into_classes(nested_paf_file)

    write_output_table(mutli_gene_dict,"multi", args.o)    
    write_output_table(single_gene_dict,"single", args.o)    

   
