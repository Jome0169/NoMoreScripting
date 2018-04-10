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
    annotation_key = []

    with open(ahrd_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split('\t')
            if "#" in clean_line[0]:
                pass
            elif clean_line == None:
                pass
            else:
                annotation_key.append(clean_line)
    return annotation_key

def read_iprscan(iprscan_file):
    """reads in interproscan output file. Creates a list for later processsing.

    :iprscan_file: interpro scan file to read
    :returns: a nested list.

    """
    ipscan_dict = {}
    with open(iprscan_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split('\t')
            
            if clean_line[0] not in ipscan_dict:
                
                if len(clean_line) > 11:
                    reformat_line =[clean_line[5],clean_line[11]]
                    ipscan_dict[clean_line[0]] = [reformat_line]
                else:
                    pass

            elif clean_line[0] in ipscan_dict:

                if len(clean_line) > 11:
                    reformat_line =[clean_line[5],clean_line[11]]
                    ipscan_dict[clean_line[0]].append(reformat_line)
                else:
                    pass

    return ipscan_dict


def read_GO_file(go_file):
    """reads in GO file, returns a list like element

    :go_file: TODO
    :returns: TODO

    """
    go_dict = {}
    with open(go_file, 'r') as f:
        for line in f:
            clean_line = line.strip().split('\t')
            if "Sequence" not in clean_line[0] and len(clean_line) > 5:
                header_name = clean_line[0]

                GO_terms = clean_line[5].split(';')
                GO_funcs = clean_line[6].split(";")
                
                #go_dict[header_name] = [first_hit]
                X = zip(GO_funcs, GO_terms)
                fixed = [list(i) for i in X]
                go_dict[header_name] = fixed

            else:
                pass
    return go_dict

def read_fasta_file(fasta_file):
    """Read in headers from fasta file. Creates a dictioninary for genes so
    that functional inforamtion can be added to the value inforamtion

    :fasta_file: fasta file
    :returns: dictionary with genes names as key and empty list as value

    """
    header_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                clean_header = line.strip().replace('>','').split()
                if clean_header[0] not in header_dict:
                    header_dict[clean_header[0]] = []
                else:
                    pass
            else:
                pass
    return header_dict

def combine_dict_AHRD(key_dict, ahrd_lines):
    """begins by appending ahrd information to the appropriate gene value in
    the gene dict

    :key_dict: TODO
    :ahrd_lines: TODO
    :returns: TODO

    """
    for item in ahrd_lines:
        if len(ahrd_lines) > 1:
            ahrd_key = item[0]
            if ahrd_key in key_dict:
                create_ahrd_list = [item[1], item[3]]
                key_dict[item[0]].append(create_ahrd_list)
            else:
                pass
        else:
            pass
    return key_dict

def combine_dict_iprscan(key_dict,ipscan_dict):
    """Appends ipscan information to the gene dict

    :key_dict: TODO
    :ipscan_dict: TODO
    :returns: TODO

    """
    for gene_name, other_info in  key_dict.items():
        if gene_name in ipscan_dict:
            key_dict[gene_name].append(ipscan_dict[gene_name])
        elif gene_name not in ipscan_dict:
            key_dict[gene_name].append([])
    
    #for ipscan_key, ipscan_info in ipscan_dict.items():
    #    if ipscan_key in key_dict:
    #        key_dict[ipscan_key].append(ipscan_info)
    #    else:
    #        key_dict[ipscan_key].append([])

    return key_dict

def combine_dict_GOterms(key_dict, GO_dict):
    """combiens the go terms found with the annotation..

    :key_dict: TODO
    :GO_dict: TODO
    :returns: TODO

    """

    for gene_name, func_info in key_dict.items():
        if gene_name in GO_dict:
            key_dict[gene_name].append(GO_dict[gene_name])
        elif gene_name not in GO_dict:
            key_dict[gene_name].append([])

        else:
            print("The following key is not found %s, program stop" % gene_name)
            sys.exit(-2)

    return key_dict 


def micro_formatting(list_reformat):
    """reforamts list of Fun Hits

    :lsit_reformat: TODO
    :returns: TODO

    """
    final_list = ""

    if len(list_reformat) > 1:
        for item in list_reformat:
            item_short = item[1]
            item_long = item[0]
            format_final = item_short + ',' + item_long + ' | '
            #add fixed string to final    
            final_list += format_final

    elif len(list_reformat) == 1:
        take_short_name = list_reformat[0][1]
        take_long_func_name = list_reformat[0][0]
        format_final = take_short_name + ',' + take_long_func_name + ' | '
        #add fixed string to final    
        final_list += format_final
    elif len(list_reformat) == 0:
        return None
    
    #We do not want a final '|' at the end 
    trim_far_bar  = final_list.rsplit('|', 1)[0] 

    return trim_far_bar 

def remove_files(file_name):
   """removes file name if it exits

   :file_name): TODO
   :returns: TODO

   """

   try:
       os.remove(file_name)
   except OSError:
       pass




def write_output_file(key_dict, output_file):
    """TODO: Docstring for write_output_file.

    :key_dict: TODO
    :output_file: TODO
    :returns: TODO

    """

    remove_files(output_file)    

    with open(output_file, 'a') as f:
        header_string = "Gene ID" + '\t' + "Human Readable Description" + '\t' \
        + 'Interpro' + '\t' "Gene Ontology" + '\n'
        f.write(header_string) 
        for gene_name, protein_info in key_dict.items():
            human_readable = protein_info[0][1]
            #Format strings for writing to the final GO term file 
            interpro_scan_string = micro_formatting(protein_info[1])
            go_term_string = micro_formatting(protein_info[2])


            
            if interpro_scan_string == None:
                interpro_scan_string = ' '
            else:
                pass

            if go_term_string == None:
                go_term_string = ' '
            elif go_term_string == ",  |":
                go_term_string = ' '
            else:
                pass
            
            final_string = gene_name + '\t' + human_readable + '\t' \
            + interpro_scan_string + '\t' + go_term_string


            #print(final_string)
            #input()
            
            f.write(final_string)
            f.write('\n')
            #print(human_readable)
            #print(test_string)
            #print(test_string2)
            #input()

            

        

def get_parser():
    parser = argparse.ArgumentParser(description='The purpose of this script it\
            to take in AHRD file files, IPRscan, as well as GO output and\
            combines them to create the final functional annoation xlxs file.')
    parser.add_argument('-a','--ahrd', help='AHRD output file', required=True, dest='a')
    parser.add_argument('-f','--fasta', help='fasta file', required=True, dest='f')
    parser.add_argument('-ip','--interpro', help='interproscan file to add to output', required=True, dest='ip')
    parser.add_argument('-go','--goterms', help='go term output file to add to AHRD output file', required=False, dest='go')
    parser.add_argument('-o','--output', help='File to Write output to. Will \
            Remove file if finds a conflicting name', required=True, dest='o')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    StartTime = datetime.now()
    args = get_parser().parse_args()
   
    #Read Files
    gene_headers = read_fasta_file(args.f)
    ahrd_file = read_ahrd_csv(args.a)
    go_term_dict = read_GO_file(args.go)
    iprscan_dict = read_iprscan(args.ip)

    #Combine Files
    combine_dict_AHRD(gene_headers,ahrd_file)
    combine_dict_iprscan(gene_headers, iprscan_dict)
    combine_dict_GOterms(gene_headers,go_term_dict)

    write_output_file(gene_headers,args.o)

    #for key, item in gene_headers.items():
    #    print(key)
    #    print(item)
    #    input()

   
   
   
   
   #digest_iprscan_values(gene_headers)


    

#Speed Things
EndTime = datetime.now()
FinalTime = EndTime - StartTime

print ("Total Time %s" % (FinalTime))

