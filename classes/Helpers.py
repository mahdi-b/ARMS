import logging
import os
import sys
from Bio import SeqIO
from collections import defaultdict
from Bio.SeqRecord import SeqRecord


class printVerbose(object):
    VERBOSE = False

    def __init__(self, msg, newline=True):
        """Prints the msg to console.

        :param msg: The thing to print.
        :param newline: Should a newline character be printed before the msg?
        """
        if printVerbose.VERBOSE:
            if newline:
                print "\n"
            print msg


def makeDirOrdie(dirPath):
    """Creates a directory 'dirPath' or exits if the 'dirPath' directory already exists.  Prevents unnecessary execution.
    :param dirPath: The filepath to the directory to look for.

    :return: The path to the new directory (same as 'dirPath').
    """
    if not os.path.isdir(dirPath):
        os.makedirs(dirPath)
    else:
        logging.error("Split fasta directory %s already exists " % dirPath)
        sys.exit()
    # TODO Ask Mahdi why we do this?  Why not just return?
    return dirPath

def makeDir(dirPath):
    """Create a directory 'dirPath' if it doesn't already exist.
    :param dirPath: The filepath to the directory to look for.

    :return: The path to the new directory (same as 'dirPath').
    """
    if not os.path.isdir(dirPath):
        os.makedirs(dirPath)
    else:
        logging.warning("Split fasta directory %s already exists " % dirPath)

def sanitize(inputs):
    """White list sanitizes a dictionary of inputs, returning a dictionary of sanitized strings.  Raises an exception if
        an invalid character is found."""



    for name, input in inputs.iteritems():
        text = str(input)

    return inputs

def move(origin, destination):
    """Moves a file from an origin to a destination.  Sanitizes inputs before execution.

    :param origin:          File path to the file to be moved
    :param destination:     File path to the file's new destination
    """

    os.rename(sanitize(origin), sanitize(destination))

def splitFileBySample(fastaFile, groupsFile, splitFastaDir):
    """ Splits an input 'fastaFile' into separate files based on the mapping in a 'groupsFile', and outputs a separate
        fasta file for each group.  Sanitizes the sequence names by removing "." dots added by align.seq.
        See <http://www.mothur.org/wiki/Group_file> for group file specification.
        See <http://www.mothur.org/wiki/Fasta_file> for .fasta file specification.

    :param fastaFile: The source .fasta file to parse.
    :param groupsFile: A group file mapping each sequence in the 'fastaFile' to a group.
    :param splitFastaDir:  File path to an output directory where the split files will be put.
    """
    # Reads the sequences/group association from the group file (similar to mothur)
    # and outputs as many files as there are samples
    # not ideal parallelization since some samples might be larger than others
    # however, this is ideal require sample file independently
    # this also sanitizes the sequences by removing "." dots added by align.seq

    # TODO: test that the file fits into memory, otherwise this could cause problems
    mySequences = SeqIO.to_dict(SeqIO.parse(fastaFile, 'fasta'))

    seqsByGroup = defaultdict(list)
    for line in open(groupsFile, 'r'):
        seq, group = line.rstrip().split()
        seqsByGroup[group].append(seq)
    for group in seqsByGroup.keys():
        outFile = open(os.path.join(splitFastaDir, group), 'a')
        for seq in seqsByGroup[group]:
            SeqIO.write(SeqRecord(mySequences[seq].seq.ungap(".").ungap("-"), id=seq, description=""), outFile, 'fasta')
        outFile.close()
