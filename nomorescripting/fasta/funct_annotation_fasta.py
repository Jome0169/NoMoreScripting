
"""
>AT1G01020.1 | ARV1 family protein | Chr1:6915-8666 REVERSE LENGTH=738 | 201606
"""

from datetime import datetime
from sys import argv
import copy
import argparse
import sys
import os

def scaf_header_load(fasta_file):
    """read in fasta file

    :fasta_file: TODO
    :returns: TODO

    """
    headers = {} 
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                clean_line = line.strip().split()
                
                gene_name = clean_line[0].replace('>','')
                other_name = ' '.join(clean_line[1:])
                headers[gene_name] = []
    
    return headers


def store_fasta_info(fasta_file):
    """TODO: Docstring for store_fasta_info.

    :fasta_file: TODO
    :returns: TODO

    """
    headers = {} 
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                clean_line = line.strip().split()
                gene_name = clean_line[0].replace('>','')
                headers[gene_name] = ''
            else:
                seqline = line.strip()
                headers[gene_name] += seqline
    return headers 

def read_in_gff(gff_file, struc, prot_file):
    """ NOTE THIS FORMAT REQUIRES THAT I USE CDS SEQ AND NOT MRNA

    :gff_file: TODO
    :returns: TODO

    """

    if struc == 'p':
        isolate_this = 'CDS'
    elif struc == 'n':
        isolate_this = 'CDS'
    elif struc == 't':
        isolate_this = 'mRNA'


    gene_info_dict = {}
    cds_all_gff_records = []
    
    with open(gff_file, 'r') as f:

        cds_seq_per_gene = []
        for line in f:
            clean_line = line.strip().split()

            if isolate_this == 'CDS':
            
                if clean_line[0].startswith("##gff-version") or \
                clean_line[0].startswith('##sequence-region'):
                    pass

                elif clean_line[0] == '###':
                    cds_all_gff_records.append(cds_seq_per_gene)
                    cds_seq_per_gene = []

                elif clean_line[2] == isolate_this:
                    cds_seq_per_gene.append(clean_line)

                else:
                    pass
            
            elif isolate_this == 'mRNA':
                if clean_line[0].startswith("##gff-version") or \
                clean_line[0].startswith('##sequence-region'):
                    pass

                elif clean_line[0] == '###':
                    cds_all_gff_records.append(cds_seq_per_gene)
                    cds_seq_per_gene = []

                elif clean_line[2] == 'exon' or clean_line[2] == \
                'three_prime_UTR' or clean_line[2] =='five_prime_UTR' or\
                clean_line[2] == 'mRNA':
                    cds_seq_per_gene.append(clean_line)

                else:
                    pass


    for gff_cds in cds_all_gff_records:
        total_seq_len = 0 
        for thing in gff_cds:
            if thing[2] == 'mRNA':
                pass
            else:
                len_of_1_cds = int(thing[4]) - int(thing[3]) 

                total_seq_len += len_of_1_cds
        
        if total_seq_len == 0:
             for thing in gff_cds:
               len_of_1_cds = int(thing[4]) - int(thing[3]) 
               total_seq_len += len_of_1_cds
        else:
            pass
        

        take_all_starts = [int(i[3]) for i in gff_cds]
        take_all_ends = [int(i[4]) for i in gff_cds]

        start_num = min(take_all_starts)
        end_num = max(take_all_ends)

        gff_cds[0][3] = start_num
        gff_cds[0][4] = end_num 
        gff_cds[0][5] = total_seq_len
        
        gene_name = gff_cds[0][8].replace('ID=','').split(';')[0].split(':')[0]
        print(gff_cds[0][8]) 
        print(gene_name)
        input()
        horrible_mRNA_len = len(prot_file[gene_name])

        gff_cds[0][5] = horrible_mRNA_len


        gene_info_dict[gene_name] = gff_cds[0]

    return gene_info_dict

    
 
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
            elif "REVERSE" in clean_line[-1] or 'FORWARD' in clean_line[-1]:
                pass
            elif clean_line == ['']:
                pass
            else:
                annotation_key[clean_line[0]] = clean_line[3]

    return annotation_key


def create_gff3_string(ID_Name,gff3_dict, struct):
    """Takes in the ID name to be looked for and the gff3_dict file with only
    mRNA headers as KEYs. Then creates a string in the TAIR genome annotaion
    format for the location info
    
    >AT1G01020.1 | ARV1 family protein | Chr1:6915-8666 REVERSE LENGTH=738 |

    :ID_Name: MRNA seq to find
    :gff3_dict: Dictionary of gff3 file with mRNA name from gff3 file as key
    :returns: Chr1:6915-8666 REVERSE LENGTH=738 

    """
   
    retrieve_gff3_annot = gff3_dict[ID_Name]

    start_loc = str(retrieve_gff3_annot[3])
    stop_loc = str(retrieve_gff3_annot[4])
    total_len = str(retrieve_gff3_annot[5])

    if retrieve_gff3_annot[6] == '+':
        strand = "FORWARD"

    elif retrieve_gff3_annot[6] == '-':
        strand = "REVERSE"


    if struct == 'p':
        prot_len = int(total_len)

        format_location_string = str(retrieve_gff3_annot[0])+':'+start_loc + '-' \
        + stop_loc + ' ' + strand + ' '+  'LENGTH='+ str(round(prot_len))

    elif struct == 'n':
        format_location_string = str(retrieve_gff3_annot[0])+':'+start_loc + '-' \
        + stop_loc + ' ' + strand + ' '+  'LENGTH='+ str(total_len)

    elif struct == 't':
        format_location_string = str(retrieve_gff3_annot[0])+':'+start_loc + '-' \
        + stop_loc + ' ' + strand + ' '+  'LENGTH='+ str(total_len)


    return format_location_string



def format_fasta_header(protein_dict, ahrd_dict, fasta_dict, gff3_dict, struc,
        location_flag):
    """Takes in all data types, and generates a final string to repalce the old
    string. 


    Wed Jun 20 14:01:19 EDT 2018: Updated this function to report locaiton
    information optinolly this assists with genomes that currently do not have
    a chromosome level assembly.
    :protein_dict: TODO
    :ahrd_dict: TODO
    :returns: TODO

    """
    final_dict = {}
    
    for gene_name, gene_info in protein_dict.items():

        if gene_name in ahrd_dict:
            functional_info = ahrd_dict[gene_name]
            
            #call to creat_gff3_string function. returns
            #Chr1:6915-8666 REVERSE LENGTH=738 
            if location_flag == True:
                #Get Location info and Len Info
                info_string = create_gff3_string(gene_name,gff3_dict, struc)
                
                functional_string = '>' + gene_name + ' ' + '|' + ' ' \
                + functional_info + ' ' + '|' + ' ' + info_string + ' ' + '|' + ' ' 

                func_String_Aed = functional_string 
                final_dict[gene_name] = func_String_Aed
            
            elif location_flag == False:
                #Create String
                functional_string = '>' + gene_name + ' ' + '|' + ' ' \
                + functional_info + ' ' + '|' + ' ' 
                
                func_String_Aed = functional_string 
                final_dict[gene_name] = func_String_Aed



    for header, seqq in fasta_dict.items():

        if header in final_dict:
            final_header = final_dict[header]
            print(final_header)
            print(seqq)
        else:
            final_header = '>' + header
            print(final_header)
            print(seqq)


def get_parser():
    parser = argparse.ArgumentParser(description='Reads in AHRD file and fasta \
            file and adds functional annotation to the FASTA. ABides by the \
            Arabidposis naming scheme')

    parser.add_argument('-f','--fasta', help='fasta file to add functional annotation', \
            required=True, dest='f')
    parser.add_argument('-ahrd','--ahrd', help='AHRD file',\
            required=True,dest='rd') 
    parser.add_argument('-g','--gff', help='gff file to read and use', \
            required=False, dest='g')
    parser.add_argument('-t','--type', help='what type of func file, optios are \
            t for trnascripts n for cds and p for protein', \
            required=True, dest='t')

    parser.add_argument('-l','--loc', help='Location added into genome name Yes|No', \
            required=False, action='store_true', dest='l')



    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    
    #ID what we're doing - protein or nuc
    prot_or_nuc = args.t 
    #Load only headers for renaming
    protein_dict_info = scaf_header_load(args.f)
    #Load Actual Fasta sequence
    fasta_seq_dict = store_fasta_info(args.f)
    
    if args.g != None:
        gff_mrna_info = read_in_gff(args.g,args.t, fasta_seq_dict)
    else:
        gff_mrna_info = None

    annotation_key = read_ahrd_csv(args.rd)
    
    format_fasta_header(protein_dict_info,annotation_key,\
            fasta_seq_dict,gff_mrna_info, args.t,args.l)



