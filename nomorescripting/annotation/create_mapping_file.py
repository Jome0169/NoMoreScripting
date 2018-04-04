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


def extract_gene_mrna_names(gff_list):
    """reads in a list of gff3 lines, takes list and parses out either mrNAs or
    the Gene names. After this creates a new name, puts them in seperate list

    :gff_list: TODO
    :returns: TODO

    """
    BaseName = "CsGy"

    mRNACounter = 10
    geneCounter = 10

    total_name_count = 6
    ChromNumber = None

    for item in gff_list:
        
        if len(item) == 2:
            pass
        else:
            if item[2] == 'mRNA':
                
                take_chrom_num = item[0].replace('Chr', '')
                if ChromNumber == None:
                    ChromNumber = take_chrom_num
                elif ChromNumber != None and take_chrom_num != ChromNumber:
                    ChromNumber = take_chrom_num
                    mRNACounter = 10
                    geneCounter = 10
                elif ChromNumber != None and take_chrom_num == ChromNumber:
                    pass

                take_chrom_num = item[0].replace('Chr', '')

                clearsplit = item[8].replace(';',':').split(':')
                real_name = clearsplit[0].replace('ID=', "")


                new_format = BaseName + str(take_chrom_num) + 'G' + str(mRNACounter).zfill(int(6)) + '.1'

                old_name_raw = real_name.split('=')
                print(real_name,'\t',new_format)
                mRNACounter += 10
            elif item[2] == 'gene':

                take_chrom_num = item[0].replace('Chr', '')
                if ChromNumber == None:
                    ChromNumber = take_chrom_num
                elif ChromNumber != None and take_chrom_num != ChromNumber:
                    ChromNumber = take_chrom_num
                    mRNACounter = 10
                    geneCounter = 10
                elif ChromNumber != None and take_chrom_num == ChromNumber:
                    pass

                take_chrom_num = item[0].replace('Chr', '')

                clearsplit = item[8].replace(';',':').split(':')
                real_name = clearsplit[1].replace('ID=', "")

                new_format = BaseName + str(take_chrom_num) + 'G' + str(mRNACounter).zfill(int(6))
                
                
                old_name_raw = real_name.split('=')
                print(old_name_raw[1],'\t',new_format)

                geneCounter += 10


read_gff = read_in_gff(argv[1])
extract_gene_mrna_names(read_gff)







