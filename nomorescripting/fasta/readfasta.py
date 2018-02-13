from datetime import datetime
from sys import argv
import copy
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys
import os
import itertools

def read_in_fasta(arg1):
    """Reads in Fasta file using biopyton

    :arg1: TODO
    :returns: TODO

    """
    file_dict = SeqIO.index(arg1, "fasta")
    return file_dict


def split_scaffold(fasta_file,scaf_name,loc_split):
    """TODO: Docstring for split_scaffold.

    :fasta_file: TODO
    :scaf_name: TODO
    :loc_split: TODO
    :returns: TODO

    """
    finallist = []
    try:
        int(loc_split)
    except ValueError:
        print("%s not a int") % loc_split

    if scaf_name in fasta_file:
        total_fasta_len = len(fasta_file[scaf_name])
        firsthalf = fasta_file[scaf_name][0:int(loc_split) - 1]
        secondhalf = fasta_file[scaf_name][int(loc_split)+1:total_fasta_len]
        
        fastastringfirst = str(scaf_name) + '_1'
        fastastringsecond = str(scaf_name) + '_2'

        finallist.append([fastastringfirst, firsthalf])
        finallist.append([fastastringsecond,secondhalf])
    
    elif scaf_name not in fasta_file:
        exit(-2)

    return finallist

def brokenscafwriter(outputname,fasta_file,scaf_name,borken_scaf,scaf_rename):
    """Writes the broken scaffold to 2 outputs. Borken scaffold, rest of
    scaffolds, 

    :outputname: basename of the broken file output
    :fasta_file: input fasta dict from biopython
    :scaf_name: scaf name that was supposed to be broken
    :borken_scaf: list of broken scaf seq

    """
    #Broken scaffold file name
    broken_scaf_file_name = str(scaf_name) +  ".broken.fasta"
    
    #If renaming string, create variable indicating if string is real
    append_to_scaf_string = None
    if scaf_name == False:
        pass
    else:
        append_to_scaf_string = scaf_rename 

    #Remove Files if they exists
    try:
        os.remove(broken_scaf_file_name)
        os.remove(outputname)
    except OSError:
        pass
    
    #At some poitn probably make a function to go through and rename
    #dict keys based off stuff. To be done another time
    with open(outputname, 'a+') as f:
        for key, val in fasta_file.items():
            if key != scaf_name and append_to_scaf_string != None:
                f.write(">" + append_to_scaf_string + str(key))
                f.write('\n')
                f.write(str(val.seq))
                f.write('\n')
            elif key != scaf_name and append_to_scaf_string == None:
                f.write(">" + str(key))
                f.write('\n')
                f.write(str(val.seq))
                f.write('\n')
            else:
                pass
    
    with open(broken_scaf_file_name, 'a+') as z:
        for item in borken_scaf:
            if append_to_scaf_string != None:
                z.write(">" + append_to_scaf_string + item[0])
                z.write('\n')
                z.write(str(item[1].seq))
                z.write('\n')
            elif append_to_scaf_string == None:
                z.write(">" + item[0])
                z.write('\n')
                z.write(str(item[1].seq))
                z.write('\n')

def rename_scaf_writer(outputname,fasta_file,scaf_rename):
    """If appending string to scaffold name, write new output fasta file. This
    seemed easier than making a bunch of if then statments below.

    :outputname: TODO
    :fasta_file: TODO
    :scaf_rename: TODO
    :returns: TODO

    """
    #If renaming string, create variable indicating if string is real
    append_to_scaf_string = None
    if scaf_name == False:                                            
        pass
    else:
        append_to_scaf_string = scaf_rename 

    #Remove Files if they exists
    try:
        os.remove(outputname)
    except OSError:
        pass

    with open(outputname, 'a+') as f:
        for key, val in fasta_file.items():
            if append_to_scaf_string != None:
                f.write(">" + append_to_scaf_string + str(key))
                f.write('\n')
                f.write(str(val.seq))
                f.write('\n')
            elif append_to_scaf_string == None:
                f.write(">" + str(key))
                f.write('\n')
                f.write(str(val.seq))
                f.write('\n')
            else:
                pass

def flaten_fasta(fasta_dict, outputname):
    """Often Times when you're playing around with other peoples Fasta, these
    things are wrapped and it makes my life a pain in the ass. So, with this I
    just wanted to write something that will simply unwrap the sequence and
    print it out.

    :fasta: TODO
    :returns: TODO

    """
    try:
        os.remove(outputname)
    except OSError:
        pass
    
    with open(outputname, 'a+') as f:
        for key, val in fasta_dict.items():
            f.write(key)
            f.write('\n')
            f.write(str(val.seq))
            f.write('\n')
 



def get_parser():
    parser = argparse.ArgumentParser(description='Software to read in fasta \
            file and do basic functionality that is often required  ')
    parser.add_argument('-f','--fasta', help='fasta file to read', \
            required=True, dest='f')
    parser.add_argument('-s','--split', help='scaffold to split',\
            required=False,dest='s') 
    parser.add_argument('-l','--locsplit', help='location in ccaffold to split',\
            required=False,dest='l') 
    parser.add_argument('-ra','--renameappend', help='append string to\
            scaffolds to rename',
            required=False,dest='ra') 
    parser.add_argument('-uw','--unwrap', help='unwrap fasta file',
            required=False,dest='uw') 

    parser.add_argument('-o','--output', help='output file to write to', \
            required=True, dest='o')
    args = vars(parser.parse_args())    
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    if args.s != None and args.l == None:
        #exit if split given without where
        print("If splitting file need location so split scaffold")
        sys.exit(-2)
    elif args.s == None and args.l != None:
        print("If splitting file need scaffold to split")
        sys.exit(-2)

    elif args.ra != None and args.s != None and args.l != None:
        file_dict = read_in_fasta(args.f)
        broken_scaffold = split_scaffold(file_dict, args.s, args.l)
        brokenscafwriter(args.o,file_dict,args.s,broken_scaffold,args.ra)
    elif args.ra != None and args.s == None and args.l == None:
        file_dict = read_in_fasta(args.f)
        rename_scaf_writer(args.o,file_dict,args.ra)
  
    #If we just want to unwrap fasta
    elif args.ra == None and args.s == None and args.l == None \
    and args.uw != None:

        file_dict = read_in_fasta(args.f)
        flaten_fasta(file_dict,args.o)


