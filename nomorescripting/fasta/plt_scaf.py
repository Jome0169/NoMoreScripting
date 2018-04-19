# -*- coding: utf-8 -*-
"""
    fasta.plt_scaf
    ~~~~~~~~~~~~~~

    plots scaffold distribution of Ns in window sizes os N. Note this is a
    module that's utalized by readfasta.py when the plt flag is used

    :copyright: (c) 2018 by YOUR_NAME.
    :license: LICENSE_NAME, see LICENSE for more details.
"""
import matplotlib.pyplot as plt
import numpy as np
from itertools import islice



def window(seq, n=2):
    "Returns a window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "

    for i in range(0, len(seq), n):
        yield seq[i:i + n] 


def calculate_chrom_len(seq_dict):
    """Calcualtes the fasta scaffold lengths.

    :seq_dict: TODO
    :returns: TODO

    """
    len_storage = {}
    for key, value in seq_dict.items():
        calc_len = len(value.seq)
        len_storage[key] = calc_len
    return len_storage



def calc_ns_in_window(seq_dict, window_size):
    """calculate the number of n's in window size of X for later use in
    plotting

    :seq_dict: TODO
    :returns: TODO

    """
    window_storage = {}
    
    if window_size == None:
        window_size = 5000
    else:
        int(window_size)
   
    for key, value in seq_dict.items():
        window_storage[key] = []
        window1 = window(value.seq, window_size)
        for item in window1:
            item = item.upper()
            #Need to edit this so that N or n 
            ncount  = (item.count('N'))
            
            window_storage[key].append(ncount)

    return window_storage


def hist_plot_ns(chrom_len_dict, ns_window_dict, window_size, scaffold_list):
    """TODO: Docstring for hist_plot_ns.

    :chom_len_dict: TODO
    :ns_window_dict: TODO
    :window_size: TODO
    :returns: TODO

    """
    if scaffold_list == None:

        for chrom, val in ns_window_dict.items():
            get_chrom_len = int(chrom_len_dict[chrom])
            create_x_range = range(0,get_chrom_len, window_size)
            #Enumerate Bins correctly
            x_pos = [i for i, _ in enumerate(create_x_range)]

            print("Plotting %s" % chrom)
            #print(len(create_x_range))
            #print(len(val))
            #print(x_pos)
            
            #Titles and width
            x_label = "Window size of %s" % str(window_size)
            y_label = "N freq, max with window size (%s)" % str(window_size)
            title = "Distribution of Ns in %s with length %s bp" % (str(chrom), \
                    str(get_chrom_len))
            bar_width = 0.5
            
            #Create and show plot
            plt.bar(x_pos,val,bar_width)
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.title(title)

            #plt.xticks(rotation=90)
            #plt.xticks(np.arange(min(x_pos), max(x_pos)+1, 1))
            creat_file_name = chrom + '_plot.pdf'
            plt.savefig(creat_file_name)
            #plt.show()
   
    #If given a list to work with, only pring and write the given ones
    elif scaffold_list != None:
        for chrom, val in ns_window_dict.items():
            if chrom in scaffold_list:
                get_chrom_len = int(chrom_len_dict[chrom])
                create_x_range = range(0,get_chrom_len, window_size)
                #Enumerate Bins correctly
                x_pos = [i for i, _ in enumerate(create_x_range)]

                print("Plotting %s" % chrom)
                #print(len(create_x_range))
                #print(len(val))
                #print(x_pos)
                
                #Titles and width
                x_label = "Window size of %s" % str(window_size)
                y_label = "N freq, max with window size (%s)" % str(window_size)
                title = "Distribution of Ns in %s with length %s bp" % (str(chrom), \
                        str(get_chrom_len))
                bar_width = 0.5
                
                #Create and show plot
                plt.bar(x_pos,val,bar_width)
                plt.xlabel(x_label)
                plt.ylabel(y_label)
                plt.title(title)

                #plt.xticks(rotation=90)
                #plt.xticks(np.arange(min(x_pos), max(x_pos)+1, 1))
                creat_file_name = chrom + '_plot.pdf'
                plt.savefig(creat_file_name)
                #Clear the plot
                plt.clf()
            else:
                pass












