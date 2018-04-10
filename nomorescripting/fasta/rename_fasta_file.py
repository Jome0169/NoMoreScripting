from sys import argv




def read_mapping_file(map_file):
    """TODO: Docstring for read_mapping_file.

    :map_file):: TODO
    :returns: dict with format of names 

    OLD - NEW
    """
    pairs = {}
    with open(map_file, 'r') as f:
        for line in f:
            clean = line.strip().split()
            pairs[clean[0]] = clean[1]
    return pairs

def GenomeReader(GenomeFile):
    """
    Arg: Takes in Genome File
    Rtrns: Returns a dictionary, Genome Scaffolds. 
    
    Keys genomic scaffold names being the
    keys - and the actual sequence being the value. 
    """
    GenomeScaffolds = {}
    with open(GenomeFile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                NamedSeq = line
                GenomeScaffolds[NamedSeq] = ""
            else:
                GenomeScaffolds[NamedSeq] += line
    return GenomeScaffolds


def rename_dict(mapping_dict, fasta_dict):
    """TODO: Docstring for rename_dict.

    :mapping_dict: TODO
    :fasta_dict: TODO
    :returns: TODO

    """

    for genename, seq in fasta_dict.items():
        isoalte_old_name = genename.split(' ')
        seq_name = isoalte_old_name[0].replace('>','')

        new_name = ">" + mapping_dict[seq_name]

        combine_new_string = new_name + ' ' + ' '.join(isoalte_old_name[1:])
        print(combine_new_string)
        print(seq)

        

read_in_genome = GenomeReader(argv[1])
map_dict = read_mapping_file(argv[2])
rename_dict(map_dict,read_in_genome)




