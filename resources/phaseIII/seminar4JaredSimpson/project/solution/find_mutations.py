from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import glob

    
# Main

# read the reference sequence
reference_filename = "EM_079517.fasta"
handle = open(reference_filename, "rU")
for record in SeqIO.parse(handle, "fasta"):
    ref_sequence = (str(record.seq))
print("ref_sequences len ", len(ref_sequence))
handle.close()


# Find all files with the name "*.genome.fasta" in the current working direectory
#  then compare each to the reference sequence and find places that they do not
#  match but ignore places where either was an 'N'

max_mutation = 0
genome_with_max = ""
min_mutation = 1000    
genome_with_min = ""
all_genome_files = glob.glob("*.genome.fasta")

mutation_dict = {}

for filename in all_genome_files:
    
    mut_count = 0
    # open filename
    handle = open(filename, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        this_sequence = str(record.seq)
    #print("this sequence len", len(this_sequence))
    # loop over and report when corresponding items are different and not 'N'
    for i in range(len(this_sequence)):
       if ((ref_sequence[i] != this_sequence[i]) and (ref_sequence[i] != 'N')\
            and this_sequence[i] != 'N'):
           print("{0} {1} {2} {3}".format(filename, i, ref_sequence[i], this_sequence[i]))
           mut_count += 1
           # put this mutation in a dictionary
           key = (i, ref_sequence[i], this_sequence[i])
           if key in mutation_dict:
               mutation_dict[key] += 1
           else:
               mutation_dict[key] = 1

    if mut_count > max_mutation:
        max_mutation = mut_count
        genome_with_max = filename
    if mut_count < min_mutation:
        min_mutation = mut_count
        genome_with_min = filename
    
    

print("genome with most mutations\n\t{0}".format(genome_with_max))
print("\t{0} mutations\n".format(max_mutation))
print("genome with fewest mutations\n\t{0}".format(genome_with_min))
print("\t{0} mutations\n".format(min_mutation))

# to find most-common mutation need to invert dictionary -- find max value and 
# corresponding key  if we don't care about ties, can just keep one mutation per count
#   a bit of a cheat

inverted = {}
for key, value in mutation_dict.items():
    inverted[value] = key

most_common_count = max(inverted.keys())
most_common_mutation = inverted[most_common_count]
print("most common mutation occurs {0} times\n".format(most_common_count))
print(most_common_mutation)


