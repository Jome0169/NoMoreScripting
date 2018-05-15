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
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde


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
    if change_consistent_cast.endswith('MB'):
        predicted_genome_size = int(change_consistent_cast.replace("MB", ''))
        final_predicted_genome_size = predicted_genome_size * 1000000

    elif change_consistent_cast.ensdwith('GB'):
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
    [N50_number, L50, over_5k, over_25k, total_assembly_size]

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

    total_assembly_size = 0
    over_5k = 0
    over_25k = 0
    seq_sizes = []
    for seq_name, val in fasta_dict.items():
        scaffold_len = (len(val.seq))
        #We only want to focus on scaffolds that appear to be real with decent
        #evidence
        if scaffold_len >= 500:
            seq_sizes.append(scaffold_len)

            if scaffold_len >= 25000:
                over_25k += 1
                over_5k += 1
            elif scaffold_len >= 5000:
                over_5k += 1

            total_assembly_size += scaffold_len
        else:
            pass
    
    N50, L50 = calculate_NX50(seq_sizes,total_assembly_size)
    #Calculate N50 using cumulative sum function in numpy
    #list_of_seq_sizes = copy.deepcopy(seq_sizes)
    #list_of_seq_sizes.sort()
    #half_of_assembly_size = total_assembly_size/2
    #cumulative_sums = np.cumsum(list_of_seq_sizes)

    ##Smallest number greater than N50 in cumulative sum list
    #find_cumulative_sum_index = min(cumulative_sums[cumulative_sums >= \
    #        half_of_assembly_size])
    #
    #index_found = np.where(cumulative_sums == find_cumulative_sum_index)
    #N50 = list_of_seq_sizes[index_found[0][0]]
    #L50 = len((cumulative_sums[cumulative_sums >= half_of_assembly_size]))
    #print(N50, L50) 

    
    #calculate #NG50 and LG50 if possiblej
    # No CLue 
    genome_prop = int(total_assembly_size) / int(genome_size)
    if genome_prop > .5:
        #Funciton call above
        NG50, LG50 = calculate_NX50(seq_sizes,genome_size)
    elif genome_prop <= .5:
        #If we can't get to 50% can't calculate N50
        NG50 = 'N/A'
        LG50 = 'N/A'
    
    final_list = [N50, L50, NG50, LG50, over_5k, over_25k, total_assembly_size]
    return final_list

def quick_read_busco_short_summary(busco_file ):
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
            elif clean_line != '':
                useful_line.append(clean_line)
    return useful_line

def plot_busco_genes(busco_list, base_name):
    """plot BUSCO genes found. I really do not wish to have to use BUSCOs janky
    R script, so i'll instead create a quick function to do it for me.

    :busco_list: busco lines read from quick_read_busco_short_summary
    :returns: TODO

    """
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

    plt.show()
    f.savefig(output_name, bbox_inches='tight')

    

def process_busco_lines(busco_list):
    """Processes busco lines and returns a list to be used in final combination
    script.

    :buso_list: TODO
    :returns: TODO

    """
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

    maker-10000558-snap-gene-0.10-mRNA-1 {'NODE_115882_length_1094_cov_3.855199': [81, 0.7901234567901234, 1.0], 
    'NODE_69614_length_1708_cov_5.135727': [81, 0.7530864197530864, 0.819672131147541]}
 

    """

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

            final_percent_ID = alignment_ID/number_hits 
            percentage_found = total_len_found/gene_total_len

            melted_scaffold_hit_dict[scaffold_name] = [gene_total_len,percentage_found,final_percent_ID]
        
        return melted_scaffold_hit_dict
            

    final_gene_dict = {}
    for gene_hit, value in nested_gene_hit_dict.items():
        melted_dict = digest_scaffold_dict(value[0])
        final_gene_dict[gene_hit] = melted_dict
    
    return final_gene_dict


def calculate_final_gene_count(melted_dict, number_CDS_genes):
    """Just counts the number of genes that pass a .90 percent ID acceptance
    threshold.

    :melted_dict: TODO
    number_CDS_genes: Number of CDS genes that were given to this file. Passed
    in by a command line argument.
    :returns: TODO


    """
    final_list = []

    gene_count = 0
    
    for genehit, scaffold_dict in melted_dict.items():
        gene_found = False
        for scaffold_name, protein_info in scaffold_dict.items():
            if float(protein_info[-1]) >= .90 and float(protein_info[1]) < 1.0:
                gene_found = True
            else:
                pass
        if gene_found == True:
            gene_count += 1 
        else:
            pass
    final_list.append(gene_count)

    if number_CDS_genes != None:
        proportian_genes_captured = int(gene_count)/int(number_CDS_genes)
        final_list.append(str(proportian_genes_captured))
    elif number_CDS_genes == None:
        final_list.append('N/A')

    return final_list


def write_output_file(base_name, assembly_stats, busco_list, gene_count_info,output_file):
    """This function writes all De-Novo information to a particular output
    file. It writes everything as a single line with the goal of all lines
    from multiple assemblies then being catted together to make analysis
    easier.

    :base_name: string of base_file_name
    :assembly_stats: list of all assembly metrics calculated by function about
    final_list = [N50_number, L50, over_5k, over_25k, total_assembly_size]

    :busco_list: TODO
    :gene_count_info: TODO
    :returns: TODO

    """
    header_string = ["base name", 'N50', 'L50','NG50',"LG50","scaffold count over 5k" , 
            "scaffold count over 25k", "total assembly size", 'Complete Busco',
            'Fragmented Busco', 'Missing Busco', "Number CDS genes found",
            "Proportion of CDS genes found"]

    #final_list = [N50_number, L50, over_5k, over_25k, total_assembly_size]


    final_string = [base_name]
    
    for item in assembly_stats:
        final_string.append(str(item))
    for count in busco_list:
        final_string.append(str(count))
    
    for values in gene_count_info:
        final_string.append(str(values))
    
    try:
        os.remove(output_file)
    except OSError:
        pass

    with open(output_file, 'a+') as f:
        #f.write(','.join(header_string))
        #f.write('\n')
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
    smalless_gene_len = []

    X = []
    Y = []
    output_name = base_name + 'genehit.scatter.pdf'

    for genehit, scaffold_dict in melted_dict.items():
        for scaffold_name, protein_info in scaffold_dict.items():
            
            
            if float(protein_info[-1]) >= .85 and float(protein_info[1]) < 1.0:
                X.append(protein_info[0])
                Y.append(protein_info[1])
                break
    
    f = plt.figure()
    #Basics Scatter
    plt.title('Gene Length Vs Percentage Recovered for %s Assembly' % base_name)
    plt.xlabel("Gene Size bp")
    plt.ylabel("Percentage of gene recovered")
    plt.scatter(X,Y)
    plt.show()
    
    #2d Histogram
    plt.hist2d(X, Y, (10, 50),cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlim(0, 10000)
    plt.xticks(range(0,10000, 1000),rotation='vertical')
    plt.show()

    ##Scatter plot with density
    xy = np.vstack([X,Y])
    z = gaussian_kde(xy)(xy)
    fig, ax = plt.subplots()
    ax.scatter(X, Y, c=z, s=100, edgecolor='')
    plt.xlim(0, 10000)
    plt.show()

    #Side by side Histogram
    create_hist_boxes = range(0,10000,500)

    #f.savefig(output_name, bbox_inches='tight')


def create_histogram_of_dups(base_name,melted_dict):
    """Many genes have multiple hits. I want to now go through and create a
    quick histogram of the number of times each gene was found. You can imagine
    that if many genes are foun 3+ times that indicative that we assembled hte
    diploid genome rather than the haploid genome.

    :melted_dict: TODO
    :returns: TODO

    """
    output_name = base_name + 'gene.duplication.hist.pdf'
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
    plt.show()
    f.savefig(output_name, bbox_inches='tight')


def create_contig_histogram(base_name,fasta_dict):
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
    
    #parser.add_argument('-head','--headers', help="Incluse headers or not in \
    #final output", required=True, dest='head')

    parser.add_argument('-genlen','--gene-number', help="Number of genes blasted \
    against in PAF file. Help calculate % genes found ", required=False, dest='gn')

    parser.add_argument('-bu','--busco', help="busco short summary output \
            file", required=False, dest='bu')
    
    parser.add_argument('-mm','--minimap2', help="Minimap 2 PAF output file to asses \
            completness", required=False, dest='bl')

    parser.add_argument('-gs','--genome-size', help="estimated genome size to \
    estimate NG50 ", required=False, dest='gs')

    parser.add_argument('-o','--output', help="output file name",\
            required=False, dest='o')
      
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    StartTime = datetime.now()
    
    ##Read Genome, will always run this.
    #fasta_file = read_in_fasta(args.f)
    #create_contig_histogram(args.base,fasta_file)
    #
    ##Calculate assemly stats. N50, and NG50 in case the genome size is included
    #if args.gs != None:
    #    genome_size = decipher_genome_size(args.gs) 
    #    assembly_stats = fasta_assembly_calculations(fasta_file,genome_size)
    #else:
    #    assembly_stats = fasta_assembly_calculations(fasta_file,None)

    #Read Busco Flad
    if args.bu != None:
        busco_lines = quick_read_busco_short_summary(args.bu)
        plot_busco_genes(busco_lines,args.base)
        final_busc_numbers = process_busco_lines(busco_lines)
    else:
        busco_lines = None
    
    #PAF_dict = parse_PAF_file(args.bl)
    #gene_scaffold_info_dict = calculate_gene_scores(PAF_dict)
    #create_histogram_of_dups(args.base,gene_scaffold_info_dict)
    #plot_genelen_percentcomp(args.base,gene_scaffold_info_dict)

    #final_gene_count = calculate_final_gene_count(gene_scaffold_info_dict,args.gn)

    #write_output_file(args.base, assembly_stats, final_busc_numbers, final_gene_count, args.o )

endtime = datetime.now()
finaltime = endtime - StartTime 

print ("Total Time %s" % (finaltime))

