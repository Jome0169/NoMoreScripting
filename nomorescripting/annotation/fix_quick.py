from sys import argv




def read_in_gff(arg1):
    """TODO: Docstring for read_in_gff.

    :arg1: TODO
    :returns: TODO

    """
    all_gff_records = []
    
    with open(arg1, 'r') as f:
        for line in f:
            if '#' in line:
                pass
            else:
                cleanline = line.strip().split()
                all_gff_records.append(cleanline)
    return all_gff_records


def read_pairs(file_name):
    """TODO: Docstring for read_pairs.

    :file_name: TODO
    :returns: TODO

    """
    pairs = []

    with open(file_name, 'r') as f:
        for line in f:
            new = line.split()
            new.append(pairs)

    return pairs

read_gff = read_in_gff(argv[1])
read_key = read_pairs(argv[2])










