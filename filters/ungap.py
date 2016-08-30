import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def ungap(file_to_clean, output_file_name, gap_char, file_type):
    """Remove the reference sequences from the MACSE files and remove the non nucleotide characters from the sequences.
       we need the database seq. names to remove them from the results files

    :param file_to_clean: Filepath to the list of sequences to clean.
    :param reference_fie: The file containing all reference sequences used to align file_to_clean.
    :param output_file_name: Filepath to where the cleaned sequences should be written.
    :return Filepath to the output file.
    """

    cleaned_seqs=[]
    output = open(output_file_name, 'w')
    i = 0
    for mySeq in SeqIO.parse(file_to_clean, file_type):
        mySeq.seq = mySeq.seq.ungap(gap_char)
        cleaned_seqs.append(mySeq)
        if i % 5000 ==0:
            SeqIO.write(cleaned_seqs, output, file_type)
            cleaned_seqs = []
        i += 1
    SeqIO.write(cleaned_seqs, output, file_type)
    return output_file_name

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Usage: file_to_clean    output_file    list_of_chars_to_remove   filetype "
    else:
        ungap(*sys.argv[1:5])
