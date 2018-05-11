# -*- coding: utf-8 -*-
"""
    assembly.take_quality_metrics
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The purpose of this script is to take in a host of assembly stats from
    various de novo assemblies and output a CSV file with basics descriptors
    that can assist you in your assment of how "good" or "bad" an assembly is.
    
    This file can take in a BUSCO set, a depth histogram from bbnorm, estimated
    genome size for N50 calculation, a fasta of contigs/scaffolds.

    :copyright: (c) 2018 by YOUR_NAME.
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
import numpy
from Bio import SeqIO
from Bio.Seq import Seq


def read_in_fasta(arg1):
    """Reads in Fasta file using biopyton

    :arg1: TODO
    :returns: TODO

    """
    file_dict = SeqIO.index(arg1, "fasta")
    return file_dict

def decipher_genome_size(genome_size):
    """genome size might come in a mb measurnment or a Gb measurnments. This
    will convert the string into something that will be useful later.

    :genome_size: size in mb or Gb
    :returns: TODO

    """

    change_consistent_cast = genome_size.upper()
    if change_consistent_cast.endwith('MB'):
        predicted_genome_size = int(change_consistent_cast.replace("MB", ''))
        final_predicted_genome_size = predicted_genome_size * 1000000

    elif change_consistent_cast.endwith('GB'):
        predicted_genome_size = int(change_consistent_cast.replace("GB", ''))
        final_predicted_genome_size = predicted_genome_size * 1000000000
    elif genome_size == None:
        final_predicted_genome_size == None
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
    [N50, NG50]

    TODO: THROW ERROR IF ASSEMBLY IS LESS THAN halg the genome size

    """
    
    total_assembly_size = 0
    over_5k = 0
    over_25k = 0
    list_of_seq_sizes = []
    for seq_name, val in fasta_dict.items():
        scaffold_len = (len(val.seq))
        #We only want to focus on scaffolds that appear to be real with decent
        #evidence
        if scaffold_len >= 500:
            list_of_seq_sizes.append(scaffold_len)
            

            if scaffold_len >= 25000:
                over_25k += 1
                over_5k += 1
            elif scaffold_len >= 5000:
                over_5k += 1

            total_assembly_size += scaffold_len
        else:
            pass
    
    #Calculate N50 using cumulative sum function in numpy
    list_of_seq_sizes.sort()
    half_of_assembly_size = total_assembly_size/2
    
    cumulative_sums = numpy.cumsum(list_of_seq_sizes)
    #Smallest number greater than N50 in cumulative sum list
    find_cumulative_sum_index = min(cumulative_sums[cumulative_sums >= \
            half_of_assembly_size])
    
    index_found = numpy.where(cumulative_sums == find_cumulative_sum_index)
    N50_number = list_of_seq_sizes[index_found[0][0]]
    L50 = index_found[0][0]
    #print(over_5k)
    #print(over_25k)
    #print(N50_number)
    #print(total_assembly_size)

    final_list = [N50_number, L50, over_5k, over_25k, total_assembly_size]
    return final_list

def quick_read_busco_short_summary(busco_file):
    """Reads in the BUSCO short summary, and reads in results for compilation
    of all results. Only reads in the short_summart.txt file

    :busco_file: TODO
    :returns: TODO

    """
    useful_line = []
    with open(busco_file,'r') as f:
        for line in f:
            clean_line = line.strip()
            if clean_line.startswith('#'):
                pass
            else:
                useful_line.append(clean_line)
    return useful_line

def parse_PAF_fiel(PAF_file):
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


    gene_hit_dict = {}
    with open(PAF_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split()
            if clean_line[0] not in gene_hit_dict:
                gene_hit_dict[clean_line[0]] = [clean_line]
            elif clean_line[0] in gene_hit_dict:
                gene_hit_dict[clean_line[0]].append(clean_line)
            else:
                pass
    
    #Function call to above
    run_this_copy = copy.deepcopy(gene_hit_dict)
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

    """

    def digest_scaffold_dict(scaffold_lvl_dict):
        """Takes in scaffold dict, and reports back the melted list with the
        alignemtn ID and total len found, etc...

        :arg1: TODO
        :returns: TODO

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
                gene_total_len += int(list1[1])
                alignment_ID += alignment_ID_list

            final_percent_ID = alignment_ID/number_hits 
            percentage_found = total_len_found/gene_total_len

            melted_scaffold_hit_dict[scaffold_name] = [gene_total_len,percentage_found,final_percent_ID]
        
        return melted_scaffold_hit_dict
            

    final_gene_dict = {}
    for gene_hit, value in nested_gene_hit_dict.items():
        melted_dict = digest_scaffold_dict(value[0])
        final_gene_dict[gene_hit] = melted_dict
    

    for key,val in final_gene_dict.items():
        print(key,val)
        input()

       
        #for item in value[2:]:

        #    gene_total_len = item[1]
        #    total_len_found = int(item[3]) - int(item[2])
        #    percentage_found = total_len_found / int(gene_total_len)


        #    alignment_ID = int(item[9]) /int(item[10])

        #    print(percentage_found,alignment_ID)


 

def get_parser():
    parser = argparse.ArgumentParser(description=" A script that takes in a\
            bunch of assembly metrics from a de novo assembly including BUSCO,\
            protein blast, and estimated genome size and reports back a series\
            of inforamtive number in csv form")
    
    #Arguments 
    parser.add_argument('-f', '--fasta', help="Genome Assembly to take stats \
            of", required=True, dest='f')

    parser.add_argument('-bu','--busco', help="Busco short summary output \
            file", required=False, dest='bu')
    
    parser.add_argument('-bl','--blast', help="Blast output file to asses \
            completness", required=False, dest='bl')

    parser.add_argument('-gs','--genome-size', help="estimated genome size to \
    estimate NG50 ", required=False, dest='gs')

    parser.add_argument('-o','--output', help="output file name",\
            required=False, dest='o')
      
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    StartTime = datetime.now()
    
    #Read Genome, will always run this.
    fasta_file = read_in_fasta(args.f)
    
    #Read and interpret genome size
    if args.gs != None:
        genome_size = decipher_genome_size(args.gs) 
        fasta_assembly_calculations(fasta_file,genome_size)
    else:
        fasta_assembly_calculations(fasta_file,None)

    #Read Busco Flad
    if args.bu != None:
        busco_lines = quick_read_busco_short_summary(args.bu)
    else:
        busco_lines = None
    
    PAF_dict = parse_PAF_fiel(args.bl)
    calculate_gene_scores(PAF_dict)

endtime = datetime.now()
finaltime = endtime - StartTime 

print ("Total Time %s" % (finaltime))

