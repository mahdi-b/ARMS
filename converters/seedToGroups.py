import sys
from Bio import SeqIO


def seedToGroups(seed_fasta, groups_file):
    """Converts a fasta file of seeds to a .groups file.  Used after clustering as part of updating groups.

    :param seed_fasta: Filepath to the fasta file containing seeds.
    :param groups_file: Filepath to the output .groups file.
    :return: Filepath to the output .groups file.
    """
    with open(groups_file, 'w') as output:
        for sequence in SeqIO.parse(open(seed_fasta,'rU'), "fasta"):
            output.write("%s\t%s\n" % (sequence.id, sequence.id))
    output.close()
    return groups_file


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: input_seed_fasta output_groups_file"
    else:
        seedToGroups(*sys.argv[1:3])