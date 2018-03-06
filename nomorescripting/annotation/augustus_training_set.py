from datetime import datetime
import subprocess
from sys import argv
import copy
import argparse
import sys
import os
from os import system
import itertools
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py
from readgff import read_in_gff2


def read_prot_file(FileToParse):
    GenomicFile = {}
    with open(FileToParse, "r") as f:
        for line in f:
            line = line.strip('\n')
            if line.startswith('>'):
                HeaderLine = line
                GenomicFile[HeaderLine] = ""
            else:
                GenomicFile[HeaderLine] += line
    return GenomicFile 


def read_in_gff2(arg1):
    """TODO: Docstring for read_in_gff.

    :arg1: TODO
    :returns: TODO

    """
    all_gff_records = []
    
    with open(arg1, 'r') as f:
        parse_smart = []
        for line in f:
            if "##gff-version" in line:
                pass
            elif "##gff-version" not in line and '#' not in line:
                cleanline = line.strip().split()
                parse_smart.append(cleanline)

            elif "##gff-version" not in line and '#' in line:
                all_gff_records.append(parse_smart)
                parse_smart = []
            elif line == '' and parse_smart != None:
                all_gff_records.append(parse_smart)
            elif line == '' and parse_smart == None:
                pass
   
    #Remove Empties
    no_empty_list = [x for x in all_gff_records if x != []]
    return no_empty_list

def read_blast_db(file_name):
    """TODO: Docstring for read_blast_db.

    :file_name: TODO
    :returns: TODO

    """
    list_of_blastdb = []
    with open(file_name, 'r') as f:
        for line in f:
            clean = line.strip()
            list_of_blastdb.append(clean)
    return list_of_blastdb

def OnlyComplete(Dict):
    
    Counter = 0
    #Remove files if they exits
    try:
        os.remove("CanidateGene.fasta")
        os.remove("CanidateList.txt")
    except OSError:
        pass
    
    #If complete Gene, Run thing
    for key, value in Dict.items():
        if "complete" in key:
            Counter += 1 
            with open("CanidateList.txt", "a+") as G:
                formatkey = str(key) +'\n'
                G.write(formatkey)

            with open("CanidateGene.fasta", 'a+') as Z:
                formatkey = str(key) +'\n'
                Z.write(formatkey)
                Z.write(value)
                Z.write("\n")
        else:
            pass
    
    print("===================================")
    print("Completed Genes from PASA Assembly were written in")
    print("CanidateGene.fasta and CanidateList.txt ")
    print("there were %s complete genes compared to %s genes total" \
            % (str(Counter), str(len(Dict.keys()))))
    print("===================================")


def create_blast_db(list_of_blastdb):
    """Creates a blastdb for every blast_db file that was given.

    :list_of_blastdb: TODO
    :returns: TODO

    """
    print("===================================")
    print("Creating Blast DB out of database files")
    print(' '.join(list_of_blastdb))
    print("===================================")
    for item in list_of_blastdb:
        subprocess.call(["makeblastdb","-in", item, "-dbtype", "prot"])

def run_blast_db(list_of_blastdb, query, outputname):
    """TODO: Docstring for run_blast_db.

    :list_of_blastdb: TODO
    :returns: TODO
    """
    print('\n')
    print("===================================")
    print("Running canidata gene fasta against blast database files")
    print("===================================")
    print('\n')
    #-outfmt "6 qseqid sseqid qlen slen length pident evalue bitscore" 
    all_new_blast_files = []
    for item in list_of_blastdb:
        if outputname == None:
            copy = item
            newstring = copy.replace('.fa', '') + '_' + "canidateGenes.blast"
            all_new_blast_files.append(newstring)
            x = subprocess.Popen(["blastp","-db", item, "-query", query, \
                "-outfmt", "6 qseqid sseqid qlen slen length pident evalue bitscore", "-out", newstring]) 
        
        elif outputname != None:
            all_new_blast_files.append(outputname)
            x = subprocess.Popen(["blastp","-db", item, "-query", query[0], \
                "-outfmt", "6 qseqid sseqid qlen slen length pident qcovs", "-out", outputname]) 
        

    x.wait()
    return all_new_blast_files


def cat_blast_files(list_blast_hit):
    """Takes a list of blast file hits and cats them all together. This new
    file will then be parsed all at once looking for individuals that hit to at
    least a single protein with 50% accuracy at least once

    :list_blast_hits: TODO
    :returns: TODO

    """
    cat_string = ' '.join(list_blast_hit)
    
    
    print("===================================")
    print("Combining all blasts %s into one file named" % cat_string)
    print("pasa_vs_databases_total.blast6")
    print("===================================")
    print('\n')

    #print(cat_string)
    system('cat %s > pasa_vs_databases_total.blast6' % cat_string)



#Functions used to pull out genes that hit to at least one other protein in the
#known databases at least 50%

def ReadBlast(Filearg):
    SequencesToKee = []

    with open(Filearg, 'r') as f:
        for line in f:
            Clean = line.strip("\n").split("\t")
            if int(Clean[3]) == int(Clean[4]) and \
            int(Clean[2]) == int(Clean[3]) and \
            int(Clean[2]) == int(Clean[4]) and \
            float(Clean[5]) > 50.00 :
                SequencesToKee.append(Clean)
    return SequencesToKee

def OnlyUniqAssmbl(BlastResult):
    UniqBlast = []
    for item in BlastResult:
        if item[0] not in UniqBlast:
            UniqBlast.append(item[0])
        else:
            pass
    return set(UniqBlast)


def WriteDictSeq(item1, item2):
    with open("StrictPossibleConsvdSeq.fasta", 'a') as f:
        f.write(item1)
        f.write('\n')
        f.write(item2)
        f.write('\n')


def WriteGoodBlast(ListofI):
    try:
        os.remove("StrictPossibleConsvdSeq.blast")
    except OSError:
        pass

    with open("StrictBlastTable_PossibleConsvdSeq.blast", 'a') as Z:
        Header = "qseqid sseqid qlen slen length pident evalue bitscore"
        SplitHead = Header.split(' ')
        JoinedHead = '\t'.join(SplitHead)
        Z.write(JoinedHead)
        Z.write('\n')
        for item in ListofI:
            Q = '\t'.join(item)
            Z.write(Q)
            Z.write('\n')


def FindAssmble(passedammbl, ProtDict):
    try:
        os.remove("StrictPossibleConsvdSeq.fasta")
    except OSError:
        pass

    for key, item in ProtDict.items():
        CleanKey = key.split(' ')
        Finished = CleanKey[0].replace('>', '')
        if Finished in passedammbl:
            WriteDictSeq(key, item)
        else:
            pass


def FilterOutSelfVSelf(FileName):
    BadBois = []
    with open(FileName, 'r') as f:
        for line in f:
            Clean = line.strip("\n").split("\t")
            if Clean[0] == Clean[1]:
                pass
            else:
                BadBois.append(Clean)
    return BadBois
    
def FilterOut80s(ListofNonSelf):
    FuckedUpIndividuals = set()
    for item in ListofNonSelf:
        TakeCovLen = TakeCovLen = (float(item[4]) / float(item[2]))

        if float(item[5]) > 70.0 and int(item[6]) > 70:
            if int(item[2]) > int(item[3]):
                FuckedUpIndividuals.add(item[1])

            elif int(item[2]) < int(item[3]):
                FuckedUpIndividuals.add(item[0])

            elif int(item[2]) == int(item[3]):
                if item[0] not in FuckedUpIndividuals and item[1] not in FuckedUpIndividuals:
                   FuckedUpIndividuals.add(item[1])
                else:
                    pass
        else:
            pass
    return FuckedUpIndividuals

def LoadProtFiles(FileToParse):
    GenomicFile = {}
    with open(FileToParse, "r") as f:
        for line in f:
            line = line.strip('\n')
            if line.startswith('>'):
                HeaderLine = line
                GenomicFile[HeaderLine] = ""
            else:
                GenomicFile[HeaderLine] += line
    return GenomicFile 


def WriteDictSeq(item1, item2):
    with open("ConsvdSeq.fasta", 'a') as f:
        f.write(item1)
        f.write('\n')
        f.write(item2)
        f.write('\n')


def FindBadAssmble(Nonpass, ProtDict):
    try:
        os.remove("ConsvdSeq.fasta")
    except OSError:
        pass

    for key, item in ProtDict.items():
        FixKey = key.replace(">",'').split(' ')
        if FixKey[0] in Nonpass:
            pass
        else:
            WriteDictSeq(key, item)





def get_parser():
    parser = argparse.ArgumentParser(description='Software to read in gff \
            file and do basic functionality \
            that is often required  ')
    parser.add_argument('-g','--gff', help='gff file to read', \
            required=False, dest='g')
    
    parser.add_argument('-f','--fasta', help='Protein fasta file to read in', \
        required=True, dest='f')
    
    parser.add_argument('-b','--blast', help='DataBases to comapre completed seqs to', \
        required=False, dest='b')

    parser.add_argument('-o','--output', help='output file to write to', \
            required=True, dest='o')
    args = vars(parser.parse_args())    
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    
    fasta_file = read_prot_file(args.f)
    OnlyComplete(fasta_file)
    blast_db_list = read_blast_db(args.b)
    create_blast_db(blast_db_list)
    all_blast_files = run_blast_db(blast_db_list,"CanidateGene.fasta", None)
    cat_blast_files(all_blast_files)


    #filter 1st blast    


    print("===================================")
    print("Filtering based off of hit %")
    print("===================================")

    kept_blast_percent = ReadBlast("pasa_vs_databases_total.blast6")
    assembly_name = OnlyUniqAssmbl(kept_blast_percent)
    FindAssmble(assembly_name,fasta_file)
    WriteGoodBlast(kept_blast_percent)


    #Blast against self
    print("===================================")
    print("BLASTING VS Self")
    print("===================================")




    output = ['StrictPossibleConsvdSeq.fasta']
    create_blast_db(output)
    slf_vs_slf = run_blast_db(output, output, "self_vs_self.blast")
        
    #create consvdSeq file
    print("===================================")
    print("Filter blast vs Self")
    print("===================================")

    input()

    selv_vs_self_filter1 = FilterOutSelfVSelf("self_vs_self.blast")
    selv_vs_self_filter2 = FilterOut80s(selv_vs_self_filter1)
    FindBadAssmble(selv_vs_self_filter2, fasta_file)



