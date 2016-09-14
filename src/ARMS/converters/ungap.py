import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def ungap(file_to_clean, output_file_name, gap_chars, file_type):
    """Removes gap characters from sequences (not sequence names) in an input fasta/fastq.

    :param file_to_clean: Filepath to the list of sequences to clean.
    :param output_file_name: Filepath indicating where to write the cleaned file
    :param gap_char_list: a list of characters to ungap
    :param file_type: Either 'fasta' or 'fastq'
    :return: The filepath to the cleaned file
    """

    cleaned_seqs=[]
    output = open(output_file_name, 'w')
    i = 0
    for mySeq in SeqIO.parse(file_to_clean, file_type):
        for gapchar in gap_chars:
            mySeq.seq = mySeq.seq.ungap(gapchar)
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
