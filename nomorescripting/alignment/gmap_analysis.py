from datetime import datetime
import copy
import argparse
import sys
import os
import itertools
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py
import logging
#Hacky as heck
sys.path.insert(0,'/Users/feilab/Programming/NoMoreScripting/nomorescripting/annotation')
from readgff import read_in_gff


def function(arg1):
    """TODO: Docstring for function.

    :arg1: TODO
    :returns: TODO

    """
    pass







def get_parser():
    parser = argparse.ArgumentParser(description="Softare that reads in gmap\
            file, and reports back some basic statistics about it. This files\
            is meant to be used with various de novo assemlbies that have CDS\
            sequenced back to them using GMAP. This way we can access the\
            overall quality of a genomic assembly ")

    parser.add_argument("-g",'--gff3', help="Gff3 file from GMAP \
to be input into file.", required=True, dest='g')

    #parser.add_argument("g",'-gff3', help="Gff3 file from GMAP \
#to be input into file.", required=True, dest='g')


    args = vars(parser.parse_args())
    return parser


if __name__ == '__main__':
    args = get_parser().parse_args()

    test_me = read_in_gff(args.g)

    for item in test_me:
        print(item)
        input()
