import operator
import sys
from Bio import SeqIO
from classes.Helpers import *
from parse.parseGroupsFileToDict import parseGroupsFileToDictOfCounts


def renameWithReplicantCounts(input_fasta, groups_file, output_fasta, filetype):
    """Covnerts a fasta and a groups file to a sorted, dereplicated, fasta named by abundance.
    Specifically, each seed in the groups file has the number of sequences it represents (including itself) appended as a
    suffix.

     e.g.

    The groups file entry:
          BALI_113_ID1  BALI_113_ID2 BALI_113_ID3

    would be named as
         >BALI_113_ID1_3  in the fasta file to show that it represents 3 sequences (itself, and two other sequences)

    :param input_fasta: Input fasta/fastq file with entries for all items in the groups file.
    :param groups_file:  Input groups file showing clustering/grouping.
    :param output_fasta: Output file path.
    :param filetype: Either 'fasta' or 'fastq'
    :return: Filepath to the output fasta.
    """

    seeds = []
    seedSizes = parseGroupsFileToDictOfCounts(groups_file)

    printVerbose("Indexing reads")
    reads = SeqIO.index(input_fasta, filetype)
    printVerbose("Done indexing reads")

    printVerbose("Renaming sequences")
    for name, count in sorted(seedSizes.items(), key=operator.itemgetter(1), reverse=True):
        s = reads[name]
        s.id = "%s_%s" % (name, count)
        s.description = ""
        seeds.append(s)
        # write in chunks
        if len(seeds) == 500000:
            SeqIO.write(seeds, open(output_fasta, 'a'), filetype)
            seeds =[]
    # write the rest of the chunk
    SeqIO.write(seeds, open(output_fasta, 'a'), filetype)
    printVerbose("Done renaming sequences")
    return output_fasta

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fasta count_file outfile"
    else:
        renameWithReplicantCounts(*sys.argv[1:4])
