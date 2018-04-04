from sys import argv







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



read_files = GenomeReader(argv[1])

for key, value in read_files.items():
    if 'mRNA-2' in key or 'mRNA-3' in key or 'mRNA-4' in key \
    or 'mRNA-5' in key or 'mRNA-6' in key or 'mRNA-7' in key \
    or 'mRNA-8' in key or 'mRNA-9' in key:
        pass
    else:
        print(key)
        print(value)
