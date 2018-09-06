from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import glob

def mutations_between_pairs(filename_a, filename_b):
    ''' (str, str) -> int
    Return the number of mutations between genomes in files named 
    genome_a and genomeb
    '''

    # open file A and read into sequence
    handle = open(filename_a, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        sequenceA = str(record.seq)
    handle.close()
    # open file B and read into sequence
    handle = open(filename_b, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        sequenceB = str(record.seq)
    handle.close()
    
    mutations = 0
    # loop over pairs of genomes and count mutations     
    for i in range(len(sequenceA)):
       if ((sequenceA[i] != sequenceB[i]) and (sequenceA[i] != 'N')\
            and sequenceB[i] != 'N'):
           mutations += 1

    return mutations


# Main


# Find all files with the name "*.genome.fasta" in the current working direectory
#  then compare each to every other

all_genome_files = glob.glob("*.genome.fasta")

for filenameA in all_genome_files:
    for filenameB in all_genome_files:
        print("{0}".format(mutations_between_pairs(filenameA, filenameB)), end=" ")
    print("")


