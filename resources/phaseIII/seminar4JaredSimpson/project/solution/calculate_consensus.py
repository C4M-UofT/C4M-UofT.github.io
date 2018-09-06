from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import glob

def get_most_frequent_base(read_sequences, pos, genome_length):
    ''' value at pos has to be one of CATGN
        create dictionary and then count frequency'''

    freq = {}
    for letter in "CATGN":
        freq[letter] = 0
    # loop over all sequences and add 
    for sequence  in read_sequences:
        nucleotide = sequence[pos]
        freq[nucleotide] += 1

    # return nucleotide for which freq[nucleotide] is largest
    max_so_far = 0
    result = "N"
    for base in "CATG":
        if freq[base] > max_so_far:
            max_so_far = freq[base]
            result = base 
    return result
    

def calculate_consensus_from_reads(input_filename):
    # Open the file containing the input reads
    handle = open(input_filename, "rU")

    # Iterate over the input reads and save their sequence in a list
    read_sequences = []
    for record in SeqIO.parse(handle, "fasta"):
        read_sequences.append(str(record.seq))

    # We initialize the genome sequence as a list of Ns (unknown nucleotide)
    genome_length = len(read_sequences[0])
    genome_sequence = ['N'] * genome_length

    # Write code in the next section to calculate the nucleotide at each position
    # of the genome by finding the most frequent base at that position in the input reads
    for i in range(genome_length):
        genome_sequence[i] = get_most_frequent_base(read_sequences, i, genome_length)


    # Write the result to disk

    # We use the name of the reads file as the root of the output file name
    out_name = input_filename.replace(".reads.fasta", ".genome.fasta")

    # Open the output file
    output_handle = open(out_name, "w")

    # Make a SeqRecord object for the consensus sequence, which we will write to the output file
    out_record = SeqRecord(Seq(''.join(genome_sequence)),
                           id=out_name,
                           description="")

    SeqIO.write(out_record, output_handle, "fasta")

#
# Main
#

# Calculate the genome sequence for a single sample
calculate_consensus_from_reads("EBOV_REDC502_MinION_GUI_Conakry_2015-07-13.reads.fasta")

# Find all files with the name "*.reads.fasta" in the current working direectory
# and run calculate_consensus_from_reads on those files
    
all_input_files = glob.glob("*.reads.fasta")
for input_filename in all_input_files:
    calculate_consensus_from_reads(input_filename)

    
