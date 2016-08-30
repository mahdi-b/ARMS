import sys
from Bio import SeqIO


def seedToNames(seed_fasta, names_file):
    """Converts a fasta file of seeds to a .names file.  Used after clustering as part of updating names.

    :param seed_fasta: Filepath to the fasta file containing seeds.
    :param names_file: Filepath to the output .names file.
    :return: Filepath to the output .names file.
    """
    with open(names_file,'w') as output:
        for sequence in SeqIO.parse(open(seed_fasta,'rU'), "fasta"):
            output.write("%s\t%s\n" % (sequence.id, sequence.id))
    output.close()
    return names_file


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: input_seed_fasta output_names_file"
    else:
        seedToNames(*sys.argv[1:3])