# -*- coding: utf-8 -*-
"""
    assembly.take_quality_metrics
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The purpose of this script is to take in a host of assembly stats from
    various de novo assemblies and output a CSV file with basics descriptors
    that can assist you in your assment of how "good" or "bad" an assembly is.
    
    This file can take in a BUSCO set, a CDS PAF file from minimap, estimated
    genome size for N50 calculation, a fasta of contigs/scaffolds.

    :copyright: (c) 2018 by PABLO MENDIETA
    :license: LICENSE_NAME, see LICENSE for more details.
"""

#import statments. Lets not reinvent the wheel
from datetime import datetime
import copy
import argparse
import sys
import os
import itertools
import statistics as s
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
import logging

def read_in_fasta(arg1):
    """Reads in Fasta file using biopyton

    :arg1: TODO
    :returns: TODO

    """

    try:
        open(arg1, 'r')
    except IOError:
        print("Error: File is not open, does not exitst")

    file_dict = SeqIO.index(arg1, "fasta")
    return file_dict

def decipher_genome_size(genome_size):
    """genome size might come in a mb measurnment or a Gb measurnments. This
    will convert the string into something that will be useful later.

    :genome_size: size in mb or Gb
    :returns: TODO

    """
    
    if genome_size == None:
        logging.info("No estimated Genome size given, Ng50 will not be calculated")
        final_predicted_genome_size = None
    else:

        change_consistent_cast = genome_size.upper()
        if change_consistent_cast.endswith('MB'):
            predicted_genome_size = int(change_consistent_cast.replace("MB", ''))
            final_predicted_genome_size = predicted_genome_size * 1000000

        elif change_consistent_cast.ensdwith('GB'):
            predicted_genome_size = int(change_consistent_cast.replace("GB", ''))
            final_predicted_genome_size = predicted_genome_size * 1000000000

    return final_predicted_genome_size

def fasta_assembly_calculations(fasta_dict, genome_size):
    """calculaes many of the base measure of hte genome assembly. Total len, #
    of contigs over 5kb and 25kb, as well as N50 and L50.
    
    Calcualtes the N50 size of the given assembly. The N50 is the size of
    the scaffold that is as the 50% mark of the genome assembly. If given a
    theoretical genome size, this function will also calculate the NG50, the
    size of the scaffold that takes you to 50% of the genome size

    :fasta_dict: Biopython type seq object
    :returns: A list of either one or two values with N50 coming first
    [N50_number, L50, over_5k, over_25k, total_assembly_size_no_500,
    raw_assembly_size]

    TODO: THROW ERROR IF ASSEMBLY IS LESS THAN halg the genome size


    """
    
    def calculate_NX50(scaffold_list, value_to_find):
        """Sub function to calculate N50 metrics, weather that be NG50 or N50.
        Since N50 is simly a metric to get to HALF of a value, calculatin it is
        relativly straigh forward and only requires a list of sizes and the
        value we're looking for.

        :scaffold_list: TODO
        :value_to_find: TODO
        :returns: TODO

        """
        take_half_size = int(value_to_find)/2
        scaffold_list.sort()

        cumulative_sums = np.cumsum(scaffold_list)
        find_cumulative_sum_index = min(cumulative_sums[cumulative_sums >= take_half_size])
    
        index_found = np.where(cumulative_sums == find_cumulative_sum_index)
        NX50_metric = scaffold_list[index_found[0][0]]
        LX50_metric = len((cumulative_sums[cumulative_sums >= take_half_size]))
            
        return (NX50_metric, LX50_metric)

    logging.info("Calculating Assembly Metrics")
    total_assembly_size = 0
    raw_assembly_size = 0
    over_5k = 0
    over_25k = 0
    seq_sizes = []
    for seq_name, val in fasta_dict.items():
        scaffold_len = (len(val.seq))
        #We only want to focus on scaffolds that appear to be real with decent
        #evidence
        raw_assembly_size += scaffold_len
        if scaffold_len >= 500:
            seq_sizes.append(scaffold_len)
            total_assembly_size += scaffold_len

            if scaffold_len >= 25000:
                over_25k += 1
                over_5k += 1
            elif scaffold_len >= 5000:
                over_5k += 1
        else:
            pass
    
    N50, L50 = calculate_NX50(seq_sizes,total_assembly_size)
        
    #calculate #NG50 and LG50 if possiblej
    if genome_size == None: #Genome size might not be specified
        logging.warning('Genome size was not vgiven and will noe be calcuted')
        genome_prop = 0
    elif genome_size != None:
        genome_prop = int(total_assembly_size) / int(genome_size)

    if genome_prop > .5:
        #Funciton call above
        NG50, LG50 = calculate_NX50(seq_sizes,genome_size)
    elif genome_prop == None:
        NG50 = 'N/A'
        LG50 = 'N/A'

    elif genome_prop <= .5:
        #If we can't get to 50% don't calculate N50
        logging.info("Cannot calculate Ng50 as assembly size != half of genome, reporting as N/A")
        NG50 = 'N/A'
        LG50 = 'N/A'

    final_list = [N50, L50, NG50, LG50, over_5k, over_25k, total_assembly_size,
            raw_assembly_size]
    
    return final_list

def quick_read_busco_short_summary(busco_file):
    """Reads in the BUSCO short summary, and reads in results for compilation
    of all results. Only reads in the short_summart.txt file

    :busco_file: TODO
    :returns: TODO

    """
    if busco_file == None:
        logging.warning("No Busco file was given, Competed and missing busco will not be calculated")

        return None

    try:
        open(busco_file, 'r')
    except IOError:
        logging.info("Busco file does not exist. Will not report")


    useful_line = []
    with open(busco_file,'r') as f:
        for line in f:
            clean_line = line.strip()
            if clean_line.startswith('#'):
                pass
            elif clean_line != '':
                useful_line.append(clean_line)
    return useful_line

def plot_busco_genes(busco_list, base_name):
    """plot BUSCO genes found. I really do not wish to have to use BUSCOs janky
    R script, so i'll instead create a quick function to do it for me.

    :busco_list: busco lines read from quick_read_busco_short_summary
    :returns: TODO

    """
    
    if busco_list == None:
        return
    else:
        pass

    X = 4
    ind = np.arange(X)
    width = 1.0
    output_name = base_name +'.busco.metrics.pdf'


    single_copy_genes = int(busco_list[2].split()[0])
    duplicated = int(busco_list[3].split()[0])
    fragmented_genes = int(busco_list[4].split()[0])
    missing_busco = int(busco_list[5].split()[0])
    
    
    f = plt.figure()
    dataset1 = np.array(single_copy_genes)
    dataset2 = np.array(duplicated)
    dataset3 = np.array(fragmented_genes)
    dataset4 = np.array(missing_busco)

    p1 = plt.bar(ind,dataset1, width, color='r')
    p2 = plt.bar(ind, dataset2, width, bottom=dataset1,color='b')
    p3 = plt.bar(ind,dataset3,width,bottom=dataset2+dataset1, color='g')
    p4= plt.bar(ind,dataset4, width, bottom=dataset1+dataset2+dataset3, color='y')
    plt.legend((p1[0], p2[0], p3, p4),
            ('Singles Copy', 'Duplicated', 'Fractured', 'Missing'))

    plt.ylabel("Number of Busco genes N=1440")
    plt.title("Busco Gene Discovery Metrics for Assembly %s" %base_name)
    
    #Turn off xticks
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    #plt.show()
    f.savefig(output_name, bbox_inches='tight')

    

def process_busco_lines(busco_list):
    """Processes busco lines and returns a list to be used in final combination
    script.

    :buso_list: TODO
    :returns: TODO

    """
    if busco_list == None:
        final_busco_list = ['N/A','N/A','N/A']
    
    elif busco_list != None:
    
        completed_genes = busco_list[1].split()[0]
        fragmented_genes = busco_list[4].split()[0]
        missing_busco = busco_list[5].split()[0]
        final_busco_list = [completed_genes,fragmented_genes,missing_busco]

    return final_busco_list


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
        logging.warning("No CDS given")
        return None
    
    logging.info("Reading in CDS info from PAf file")
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





def read_gmap_PSL(PSL_file):
    """Reads in PSL file. Reads this in to a list object, and returns this list
    for later processing. 

    :PSL_file: TODO
    :returns: TODO

    """
    
    PSL_list = []
    try:
        open(PSL_file, 'r')
    except:
        logging.warning("No CDS file given")
        return None

    with open(PSL_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split()
            PSL_list.append(clean_line)
    return PSL_list

def create_gmap_gene_dict(PSL_list):
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

def gmap_nest_target_gene_dict(query_dict, PSL_list):
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


def gmap_find_split_genes(query_dict):
    """Takes in query dict and looks for geneA that might be split in between
    two different scaffolds. If this function finds evidence for this, it will
    generate a false string taht has the proper coordinates.
    
    gala12g28350    339     2       182
    gala12g28350    339     182     339

    :query_dict: TODO
    :returns: scaffold1___scaffold2 [NOne, none, none,2,339]

    """
    x_dict = copy.deepcopy(query_dict)
    split_gene_dict = {}
    for key,val in x_dict.items():
        
        get_key = val.keys()
        if len(get_key) == 2:
            #Creates iterator pairs for easy regtrieval with lambda
            x = list((itertools.combinations(get_key,2)))
            
            #Lambda function. Reports back query hit info
            query_hits = lambda dict1, key1: dict1[key1][0][10:13]

            gene_hit1 = (query_hits(val,x[0][0]))
            gene_hit2 = (query_hits(val,x[0][1]))

            if gene_hit1[2] == gene_hit2[1]:

                gene_mismatch = lambda dict1, key1: dict1[key1][0][1]
                
                gene1_mismatch = gene_mismatch(val,x[0][0])
                gene2_mismatch = gene_mismatch(val,x[0][1])

                final_scaffold_name = x[0][0] + '____' + x[0][1]
                final_query_len = gene_hit1[0]
                final_query_start = gene_hit1[1]
                final_query_end = gene_hit2[2]
                final_gene_mismatch = (int(gene1_mismatch[0]) + int(gene2_mismatch[0]))
                
                fake_list = [None] * 15

                fake_list[1] = str(final_gene_mismatch)
                fake_list[10] = str(final_query_len)
                fake_list[11] = str(final_query_start)
                fake_list[12] = str(final_query_end)

                mini_dict_scaffold = {final_scaffold_name:[fake_list]}
                
                split_gene_dict[key] = mini_dict_scaffold

            else:
                split_gene_dict[key] = val

                
        else:
            #Genes with 2+ hits are too messy and don't work well, and we don't
            #care about genes with just a single hit
            split_gene_dict[key] = val
    
    return split_gene_dict


def gmap_clobber_inner_list(neseted_query_dict):
    """Takes in nested PSL dictionary, and clobbers the most inner list into a
    string that can then actualyy be used to graph the overall amount of CDS
    retrieved. This has two advantages. Namely we don't have to parse this ugly
    ass innter list later, and it will make calculating overall CDS quality
    much easier.

    input_inner_list: ['2086', '5', '0', '0', '1', '6', '1', '264', '-', 'gala03g06192', '2118', '0', '2097', 'NODE_29262_length_3826_cov_23.417445', '3826', '1236', '3591', '3', '33,1760,298,', '21,54,1820,', '1236,1533,3293,']]
    :neseted_: {query: {target1:[hit1, hit2, hit3], target2:[hit1, hit2, hit3], etc...}}

    :returns: Something nice

    """

    def gmap_clobber_inner_inner(nested_list):
        """Sometimes a single gene will be split on the same scaffold due to
        h

        intronic sequence. this inner function takes a list, and if the lsit
        len > 1, it takes the average of PID and reports back an updated
        alignment length. 

        :nested_list): [[length_of_query,alignment_length, percent_ID,],
        [length_of_query,alignment_length, percent_ID, ]]


        :returns: TODO
        single: [length_of_query,alignment_length, percent_ID, ]


        """
        if len(nested_list) >= 2:
            gene_len = list(map(lambda x: int(x[0]),nested_list))
            avg_gene_len = sum(gene_len)/len(gene_len)

            algn_vals = list(map(lambda x: int(x[1]),nested_list))
            avg_algn_vals = sum(algn_vals)
            
            PID_vals = list(map(lambda x: float(x[2]),nested_list))
            average_PID = float(sum(PID_vals))
            
            final_clobbered_list = [avg_gene_len, avg_algn_vals,average_PID]
            nested_list = final_clobbered_list

        elif len(nested_list) == 1:
            remove_nesting = nested_list[0]
            nested_list = remove_nesting
        return nested_list 

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
            
            #Check for gene split in same scaffold
            sam_scaf_multi_hit  = gmap_clobber_inner_inner(reformulat_nested_list) 
            nested_dict[scaffold] = sam_scaf_multi_hit 

    return copy_dict

def seperate_into_classes(melted_dict):
    """Takes in melted dict, and created two seperate list based off how many
    scaffolds each gene hits to. I am curious to see how many genes hit
    multiple times, as well as how many genes hit once. Should be rather
    interesting.

    :melted_dict: TODO
    :returns: TODO

    """
    logging.info("Identifying CDS sequences found once, or more than once")
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

def calculate_final_gene_count(single_genes, multi_genes, number_CDS_genes):
    """Takes in two different dictionaries, one with sinlge genes hitting and
    one with genes that were found multiple times. Counts number of genes that
    pass a certain threshhold, and adds them to s certain value. 
    Finally counts the total proportion of number_CDS_genes captured in total.

    :melted_dict: TODO
    number_CDS_genes: Number of CDS genes that were given to this file. Passed
    in by a command line argument.
    :returns: Final_singls_cds_count, final_multi_cds_count, total_number_genes

    """
    logging.info("Calculating final CDS statistics") 
    #Overwrite list if values non empty
    final_list = ['N/A', 'N/A', 'N/A']
    single_gene_count = 0
    multi_gene_count = 0

    if single_genes == None and multi_genes == None:
        return final_list

    else:

        for gene, hits in single_genes.items():
            for scaffold_nam, hit_info in hits.items():
                if float(hit_info[-1]) > .8:
                    single_gene_count += 1 
                    break
                else:
                    pass
        
        for gene2, hits2 in multi_genes.items():
            for scaffold_nam, hit_info in hits2.items():
                if float(hit_info[-1]) > .8:
                    multi_gene_count += 1 
                    break
                else:
                    pass
    
    if number_CDS_genes != None:
        total_genes_found = int(single_gene_count) + int(multi_gene_count)
        prop_found = (total_genes_found/int(number_CDS_genes))
        final_list = [str(single_gene_count), str(multi_gene_count), str(prop_found)]

    else:
        final_list = [str(single_gene_count), str(multi_gene_count), 'N/A']

    return final_list


def write_output_file(base_name, assembly_stats, busco_list, gene_count_info,
        head_flag, output_file):
    """This function writes all De-Novo information to a particular output
    file. It writes everything as a single line with the goal of all lines
    from multiple assemblies then being catted together to make analysis
    easier.

    :base_name: string of base_file_name
    :assembly_stats: list of all assembly metrics calculated by function about
    final_list = [N50_number, L50, over_5k, over_25k, total_assembly_size,
    raw_assembly_size]

    :busco_list: TODO
    :gene_count_info: TODO
    :returns: TODO

    """

    header_string = ["base name", 'N50', 'L50','NG50',"LG50","scaffold count over 5k" , 
            "scaffold count over 25k", "total assembly size scaffolds over 500", 
            'Raw Assembly size', 'Complete Busco',
            'Fragmented Busco', 'Missing Busco', 'Number of Single CDS found', 
            'Number of CDS found once', "Proportion of CDS found"]
    #final_list = [N50_number, L50, over_5k, over_25k, total_assembly_size]


    final_string = [base_name]
    
    for item in assembly_stats:
        final_string.append(str(item))
    for count in busco_list:
        final_string.append(str(count))
    
    for values in gene_count_info:
        final_string.append(str(values))
    
    try:
        #Should be a log message here.
        logging.warning("Output file %s will be removed." % output_file)
        os.remove(output_file)
    except OSError:
        pass

    with open(output_file, 'a+') as f:
        if head_flag == True:
            f.write(','.join(header_string))
            f.write('\n')
        else:
            print("No Header Inforation will be added to output")
            pass
        f.write(','.join(final_string))
        f.write('\n')


def plot_genelen_percentcomp(base_name,melted_dict):
    """Some of the genes retrieved were not of the highest quality. This makes
    sense just due to the nature of gene alignment. This script goes through
    and narrows down a list based off of a relative gene %ID of 60% or greater.
    I will then go through and create a scatter plot of all these genes to
    asses their overall quality.

    :melted_dict: most inner layer has the format of 
    melted_scaffold_hit_dict[scaffold_name] = [gene_total_len,percentage_found,final_percent_id]
        
    :returns: TODO

    """
    if melted_dict == None:
        return
    else:
        pass

    smalless_gene_len = []

    X = []
    Y = []
    output_name = base_name + '.genehit.scatter.pdf'

    for genehit, scaffold_dict in melted_dict.items():
        for scaffold_name, protein_info in scaffold_dict.items():
            
            
            if float(protein_info[-1]) >= .85 and float(protein_info[1]) < 1.0:
                X.append(protein_info[0])
                Y.append(protein_info[1])
                break
    
    #f = plt.figure()
    ##Basics Scatter
    #plt.title('Gene Length Vs Percentage Recovered for %s Assembly' % base_name)
    #plt.xlabel("Gene Size bp")
    #plt.ylabel("Percentage of gene recovered")
    #plt.scatter(X,Y)
    ##plt.show()
    #f.savefig(output_name, bbox_inches='tight')

    ##2d Histogram
    #plt.hist2d(X, Y, (10, 50),cmap=plt.cm.jet)
    #plt.colorbar()
    #plt.xlim(0, 10000)
    #plt.xticks(range(0,10000, 1000),rotation='vertical')
    ##plt.show()
    #f.savefig(output_name, bbox_inches='tight')

    ##Scatter plot with density
    xy = np.vstack([X,Y])
    z = gaussian_kde(xy)(xy)
    fig, ax = plt.subplots()
    ax.scatter(X, Y, c=z, s=100, edgecolor='')
    plt.xlim(0, 10000)
    plt.title('Gene Length Vs Percentage Recovered for %s Assembly' % base_name)
    plt.xlabel("Gene Size bp")
    plt.ylabel("Percentage of gene recovered")

    plt.draw()
    #plt.show()

    fig.savefig(output_name, bbox_inches='tight')

def plot_histogram_of_dups(base_name,melted_dict):
    """Many genes have multiple hits. I want to now go through and create a
    quick histogram of the number of times each gene was found. You can imagine
    that if many genes are foun 3+ times that indicative that we assembled hte
    diploid genome rather than the haploid genome.

    :melted_dict: TODO
    :returns: TODO

    """
    if melted_dict == None:
        return
    else:
        pass

    output_name = base_name + '.gene.duplication.hist.pdf'
    count_list = []

    for key, value in melted_dict.items():
        count_list.append(len(value))
    
    f = plt.figure()
    bins = np.arange(min(count_list),20,1.0)
    plt.xlim(0, 20)
    plt.xticks(range(0,20))
    
    plt.title('Frquency of finding gene X times in %s Assembly' % base_name)
    plt.xlabel("Number of Times Gene Found")
    plt.ylabel("Frequency of Genes being Found")

    plt.hist(count_list,bins=bins,alpha=.5,edgecolor='black', linewidth=1.0)
    #plt.show()
    f.savefig(output_name, bbox_inches='tight')


def plot_contig_histogram(base_name,fasta_dict):
    """Creates a histogram of all the scaffolds and their sizes with bin sizes
    scaled. This just works as a great graphical metric to get an idea of how
    the assembly is behaving.

    :fasta_dict: TODO
    :returns: TODO

    """
    output_name = base_name + '.contig.hist.pdf'
    contig_len_list = []
    for seq_name, val in fasta_dict.items():
        contig_len = (len(val.seq))
        if contig_len > 500:
            contig_len_list.append(contig_len)
        else:
            pass

    f = plt.figure()
    bins = np.arange(500,max(contig_len_list),500)

    plt.xticks(range(0,max(contig_len_list),1000), rotation='vertical')
    plt.title('histogram of contig lengths for %s genome' % base_name)
    plt.xlabel("Contig Size")
    plt.ylabel("Frequency of Contig Size")
    plt.hist(contig_len_list,bins=bins,alpha=.5,edgecolor='black', linewidth=.3)
    #plt.show()

    f.savefig(output_name, bbox_inches='tight')
 
 

def get_parser():
    parser = argparse.ArgumentParser(description=" A script that takes in a\
            bunch of assembly metrics from a de novo assembly including BUSCO,\
            protein blast, and estimated genome size and reports back a series\
            of inforamtive number in csv form")
    
    #Arguments 
    parser.add_argument('-f', '--fasta', help="Genome Assembly to take stats \
            of", required=True, dest='f')

    parser.add_argument('-base','--base-name', help="Base name of the file to \
    be used to the first cell of the final output.", required=True, dest='base')
    
    parser.add_argument('-head','--headers', help="Incluse headers or not in \
    final output", action='store_true', required=False, dest='head')

    parser.add_argument('-genelen','--gene-number', help="Number of genes blasted \
    against in PAF file. Help calculate percentage genes found ", required=False, dest='gn')

    parser.add_argument('-bu','--busco', help="busco short summary output \
            file", required=False, dest='bu')
    
    parser.add_argument('-mm','--minimap2', help="Minimap 2 PAF output file to asses \
            completness", required=False, dest='mm')

    parser.add_argument('-gmap','--gmap-psl', help="gmap file output for cds seq \
            completness", required=False, dest='gmap')


    parser.add_argument('-gs','--genome-size', help="estimated genome size to \
    estimate NG50 ", required=False, dest='gs')

    parser.add_argument('-o','--output', help="Output file name",\
            required=False, dest='o')
      
    return parser

if __name__ == "__main__":
    
    args = get_parser().parse_args()
    StartTime = datetime.now()

    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p:')
    logging.info("Starting %s script" % __file__)

    #Read Genome, will always run this.
    logging.info("Reading in Fasta file: {}".format(args.f))
    fasta_file = read_in_fasta(args.f)
    
    logging.info("Plotting kmer histogram for: {}".format(args.f))
    plot_contig_histogram(args.base,fasta_file)
    
    #Calculates genome size based off ending of value mb or gb.
    #returns either value or None is not included
    genome_size = decipher_genome_size(args.gs) 

    #Calculate assemly stats. N50, and NG50 in case the genome size != None
    assembly_stats = fasta_assembly_calculations(fasta_file,genome_size)

    #Read Busco file
    busco_lines = quick_read_busco_short_summary(args.bu)
    plot_busco_genes(busco_lines,args.base)
    final_busc_numbers = process_busco_lines(busco_lines)
       

    #Read in CDS sequence. Based on if minimap, or gmap was used
    
    if args.mm != None and args.gmap == None:
        logging.info("Minimap2 output PAF file included. Calculating CDS")
        #Read in PAF file in dict format
        PAF_dict = parse_PAF_file(args.mm)

        #Nest Genes
        nested_paf_file = calculate_gene_scores(PAF_dict)

        #function in this file
        single_gene_dict,mutli_gene_dict = seperate_into_classes(nested_paf_file)
        
        #Plot found genes
        #plot_histogram_of_dups(args.base,gene_scaffold_info_dict)
        #plot_genelen_percentcomp(args.base,gene_scaffold_info_dict)
        
        #caculate final # of genes found with percentage cut off >.6
        final_gene_count = calculate_final_gene_count(nested_paf_file,args.gn)

    elif args.mm == None and args.gmap != None:
        #Read in gmap file
        logging.info("GMAP output PSL file included. Calculating CDS capture")

        gmap_list = read_gmap_PSL(args.gmap)
        #create Dict
        gmap_dict = create_gmap_gene_dict(gmap_list)
        #Nest by scaffold hit
        gmap_nest_target_gene_dict(gmap_dict, gmap_list)
        #Find split genes
        gmap_fixed_split_genes = gmap_find_split_genes(gmap_dict )
        #Melt for PID and gene len capture
        gmap_melted_dict = gmap_clobber_inner_list(gmap_fixed_split_genes)

        single_gene_dict,mutli_gene_dict = seperate_into_classes(gmap_melted_dict)

        #Plot found genes
        #plot_histogram_of_dups(args.base,gene_scaffold_info_dict)
        #plot_genelen_percentcomp(args.base,gene_scaffold_info_dict)



        #plot_histogram_of_dups(args.base,gmap_melted_dict)
        final_gene_count = calculate_final_gene_count(single_gene_dict, \
                mutli_gene_dict,args.gn)

    elif args.mm == None and args.gmap == None:

        logging.info("No CDS alignment file included. CDS capture will not be \
                reported")
        final_gene_count = calculate_final_gene_count(None,None)

    #Write output
    write_output_file(args.base, assembly_stats, final_busc_numbers,
            final_gene_count, args.head ,args.o )

    logging.info("Script done! Find output in %s" % args.o)

endtime = datetime.now()
finaltime = endtime - StartTime 
logging.info("Total length of time for run was %s" % finaltime)

