import getopt
import os
from sys import argv

"""
Takes in the completed files from pasa training set creator and only outputs
individuals with "COMPLETE" in their name from the PEP file. These individuals
have their sequence written to canidatelist.txt and canidatelist.fasta

"""





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




def OnlyComplete(Dict):
    for key, value in Dict.iteritems():
        if "complete" in key:
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

try:
    os.remove("CanidateGene.fasta")
    os.remove("CanidateList.txt")
except OSError:
    pass



LoadedFile = LoadProtFiles(argv[1])
OnlyComplete(LoadedFile)







#def Usage():
#
#def Main():
#
#
#    try:
#        options, other_args = getopt.getopt(sys.argv[1:], "i:b:h:o:", ["help"])
#
#
#
#
#    except getopt.GetoptError:
#        print "There was a command parsing error"
#        Usage()
#        sys.exit(1)








