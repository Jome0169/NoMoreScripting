# -*- coding: utf-8 -*-
"""
    assembly.PAF_reader
    ~~~~~~~~~~~~~~~~~~~

    Basic file to read in PAF alignment format from minimap2. Minimap2 is much
    faster and better than blast, so might as well figure out a good way to
    read and write such information. I would eventually like this script to be
    able to be imported into whatever script that can input PAF file format.

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
    
    #Below it the filter low mappign quality command. THis has been deprecated
    #as i've realize that filtereing by mapping quality doesn't make send when
    #looking for presence abscene data
    #def filter_low_mapping_qual(scaffold_lvl_dict):
    #    """Minimap2 asssigns a mapping quality metric when it creates
    #    algignemtns. This can be useful when trying to identify canidate CDS
    #    geens that aligned anc actually perform well. Therefore, this functiona
    #    aims to clear out any alignemtsn with less than a .999 certainy in
    #    alignment.

    #    :scaffold_lvl_dict: TODO
    #    :returns: TODO

    #    """
    #    edited_scaffold_lvl_dict = copy.deepcopy(scaffold_lvl_dict)
    #    cleaned_dict_no_low_map = {}

    #    for scaffold_name, hits in edited_scaffold_lvl_dict.items():
    #        good_hits = []
    #        for list1 in hits:
    #            if int(list1[11]) >= 30:
    #                good_hits.append(list1)
    #            elif int(list1[11]) < 30:
    #                pass
    #        cleaned_dict_no_low_map[scaffold_name] = good_hits

    #    return cleaned_dict_no_low_map


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
            #high_qual_mapping_dict = filter_low_mapping_qual(value[0])
            melted_dict = digest_scaffold_dict(value[0])
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


