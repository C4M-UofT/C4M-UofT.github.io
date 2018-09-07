from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

filename = "EBOV_REDC502_MinION_GUI_Conakry_2015-07-13.reads.fasta"

# Open the file containing the input reads
handle = open(filename, "rU")

# Iterate over the input reads and save their sequence in a list
read_sequences = []
for record in SeqIO.parse(handle, "fasta"):
    read_sequences.append(str(record.seq))

print(len(read_sequences))
print(len(read_sequences[0]))
