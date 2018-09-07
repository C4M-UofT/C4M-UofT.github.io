from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
import glob

output_handle = open("all_genomes.fasta", "w")
for filename in glob.glob("*.genome.fasta"):
    input_handle = open(filename, "rU")
    record = SeqIO.parse(input_handle, "fasta")
    SeqIO.write(record, output_handle, "fasta")
output_handle.close()
