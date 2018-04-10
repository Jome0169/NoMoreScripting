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

def read_in_gff(gff_file):
    """TODO: Docstring for read_in_gff.

    :gff_file: TODO
    :returns: TODO

    """
    all_gff_records = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            cleanline = line.strip().split('\t')
            all_gff_records.append(cleanline)
    return all_gff_records



def fix_gene(gene_string, map_file_dict):
    """todo: docstring for fix_gene.
    id=csgy0g000020;name=csgy0g000020;note=unknown protein

    """
    
    final_string = ""
    
    split_gene = gene_string.split(';')

    rejoin_w_colon = []
    for item in split_gene:
        split_on_equals = item.split('=')
        if split_on_equals[1] in map_file_dict:
            find_new_gene = map_file_dict[split_on_equals[1]]
            split_on_equals[1] = find_new_gene
        else:
            pass
        
        rejoin_w_colon.append('='.join(split_on_equals))

    final_string = ';'.join(rejoin_w_colon)

    return final_string



def fix_other(gene_string, map_file_dict):
    """todo: docstring for fix_gene.
    
    ID=CsGy0G000030.1:cds;Parent=CsGy0G000030.1
    ID=CsGy0G000030.1:three_prime_utr;Parent=CsGy0G000030.1
    """
    
    
    split_gene = gene_string.split(';')
    rejoin_w_semicolon = []
    for item in split_gene:
        rejoin_w_equals = [] 
        split_on_equals = item.split('=')
        for item in split_on_equals:
            split_colons  = []
            if ':' in item:
                once_more = item.split(':')

                #if once_more[0] in map_file_dict:
                #    correct_name = mapping_dict[once_more[0]]
                #    split_colons.append(correct_name)
                #else:
                #    split_colons.append([1:])

                for thing1 in once_more:
                    if thing1 in map_file_dict:
                        correct_name = map_file_dict[thing1]
                        split_colons.append(correct_name)
                    else:
                        split_colons.append(thing1)
                rejoin_w_real_colon = ':'.join(split_colons) 
                rejoin_w_equals.append(rejoin_w_real_colon)
            
            elif ':' not in item:
                if item in map_file_dict:
                    new_name = map_file_dict[item]
                    rejoin_w_equals.append(new_name)
                else:
                    rejoin_w_equals.append(item)
        rejoin_w_semicolon.append('='.join(rejoin_w_equals))
                    
    final_string =(';'.join(rejoin_w_semicolon))
    return final_string

def fix_function(lines_gff3, map_file_dict):
    """TODO: Docstring for function.

    :arg1: TODO
    :returns: TODO

    """
    for line in lines_gff3:
        if line[0].startswith('#'):
            print(' '.join(line))
        
        elif line[2] == 'gene':
            fixed_gene_line = fix_gene(line[8], map_file_dict)
            line[8] = fixed_gene_line 
            print('\t'.join(line))

        elif line[2] == 'mRNA':
            fixed_mrna_line = fix_gene(line[8], map_file_dict)
            line[8] = fixed_mrna_line
            print('\t'.join(line))
        else:
            fixed_other = fix_other(line[8], map_dict)
            line[8] = fixed_other
            print('\t'.join(line))






gff_list = read_in_gff(argv[1])
map_dict = read_mapping_file(argv[2])
fix_function(gff_list,map_dict)







