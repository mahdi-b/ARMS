from Bio import SeqIO
import sys
import os

"""
Takes in a fasta file and outputs a new fasta with the sequences renamed.  Renaming convention for sequences is
x<n> for x.y.z.fasta, where n is an integer in the range [0:n] where n is the position of the sequence in the
input_file.

e.g.
# generates a_renamed.fasta, with sequence names as aID0, aID1, aID2 ...
python rename_sequences.py a.b.c.fasta a_renamed.fasta fasta

rename("a.b.c.fasta, "a_renamed.fasta", "fasta")

"""

def serialRename(input_file, output_file, file_type):
    """Takes in a fasta file and outputs a new fasta with the sequences renamed.  Renaming convention is x.y.z<n> for
        x.y.z.fasta, where n is an integer in the range [0:n] where n is the position of the sequence in the input_file.

    :param input_file:      Input fasta or fastq file.
    :param output_file:     Filepath for the output file.
    :param file_type:       "fasta" or "fastq"
    :return:
    """
    prefix = os.path.splitext(input_file)
    if input_file.split('.')[0]:
        i = 0
        with open(input_file, 'r') as input:
            with open(output_file, 'w') as output:
                for s in SeqIO.parse(input, file_type):
                    s.id ="%sID%s" % (prefix, i)
                    s.description = ""
                    SeqIO.write(s, output, file_type)
                    i += 1


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: input_file output_file renaming_prefix file_type"
        exit()
        serialRename(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print "Could not open input_file."
        exit()