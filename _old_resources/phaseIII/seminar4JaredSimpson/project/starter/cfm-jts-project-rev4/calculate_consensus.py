from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import glob

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

    # Write code in the next section to calculate the nucleotide at each 
    # position of the genome by finding the most-frequent base at that position
    # in the input reads
    # HINT: You will probably find it helpful to create a helper function that 
    # does the task of finding the most-frequent base for a particular position.
    # Then call your function for each position. 

    #
    # YOUR CODE STARTS HERE
    #

    #
    # YOUR CODE ENDS HERE
    #

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
    
#
# YOUR CODE HERE
#
