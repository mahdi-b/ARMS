import logging
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from multiprocessing import Pool

class printVerbose(object):
    """A class to toggle verbose printing."""

    # If True, verbose printing is enabled
    VERBOSE = False

    def __init__(self, msg, newline=True):
        """Prints the msg to console if printVerbose.VERBOSE is True.
        :param msg: The thing to print.
        :param newline: Should a newline character be printed before the msg?
        """
        if printVerbose.VERBOSE:
            if newline:
                print "\n"
            print msg

class MissingMothurFilterException(Exception):
    def __str__(self):
        return "Expected a mothur filter parameter, but found none"

class MissingMothurFileException(Exception):
    def __str__(self):
        return "Expected a mothur file to update, but found none."


def runInstance(args, myInstance):
    """Runs an instance of a ProgramRunner.  Calls ProgramRunner.run() for a ProgramRunner object.'
        :param args:        Arguments from the argparse object.
        :param myInstance A fully initalized ProgramRunner object to run.
    """
    # Use the global version to facilitate calling workers
    print "runinstance"
    if args.dryRun:
        logging.info(myInstance.dryRun())
    else:
        # logging.info(myInstance.dryRun())
        myInstance.run()




def parallel(function, args, runners, pool=Pool(processes=1)):
    '''
    Executes one or more ProgramRunners in parallel.
    :param function:    The function to call over ProgramRunners.  Generally runInstance().
    :param args:        Arguments from the argparse object.
    :param runners:     A list of fully initalized ProgramRunners to execute.
    :param pool:        An initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    :return:
    '''

    data = [(args, runner) for runner in runners]
    for arg, runner in data:
        runInstance(arg, runner)
        # pool.map(function, data)

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
    # TODO Ask Mahdi why we do this?  Why not just return?  See makeDir below...
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
        an invalid character is found.
    :param inputs: The list of inputs to sanitize

    :return: A list of sanitized inputs
    """

    for name, input in inputs.iteritems():
        text = str(input)

    return inputs



def move(origin, destination):
    """Moves a file from an origin to a destination.  Sanitizes inputs before execution.

    :param origin:          File path to the file to be moved.  e.g. "dirA/a.txt"
    :param destination:     File path to the file's new destination.  e.g. "dirB/a.txt"
    """
    os.rename(sanitize(origin), sanitize(destination))

"""
def mothur_parseInputFileType(args):
    Gets the input file format and input file from the args argparser.

    :param args:
    :return: A tuple: (input file format, input file path)

    if args.fasta:
        inputFileType = "fasta"
        inputFile = args.fasta
    elif args.list:
        inputFileType = "list"
        inputFile = args.list
    elif args.groups:
        inputFileType = "groups"
        inputFile = args.groups
    elif args.names:
        inputFileType = "names"
        inputFile = args.names
    elif args.count:
        inputFileType = "count"
        inputFile = args.count

    #TODO not sure if these are the right keywords
    elif args.contigsReport:
        inputFileType = "contigsreport"
        inputFile = args.contigsReport

    elif args.summaryFile:
        inputFileType = "summary"
        inputFile=args.summaryFile
    else:
        # alignReport
        inputFileType = "alignreport"
        inputFile = args.alnReport

    return inputFileType, inputFile
"""


def mothur_buildOptionString(args, mustFilter=False, mustUpdate=False):
    """Builds a mothur parameter string (the optional filter/update section of a mothur command) form the argparser \
        using a known list of filters and filetypes.

    :param args: An argparse object with...
        any of the following filter parameters:
        start	        Maximum allowable sequence starting index.
        end	            Minimum allowable sequence ending index.
        minlength	    Minimum allowable sequence length.
        maxlength	    Maximum allowable sequence length.
        maxambig	    Maxmimum number of allowed ambiguities.
        maxn	        Maximum number of allowed N's.
        maxhomop	    Maximum allowable homopolymer length.
        and any of the following update files...

        groups  	    Groups file to be updated.  See <http://www.mothur.org/wiki/Group_file>
        names           Names file to be updated.  See <http://www.mothur.org/wiki/Name_file>
        alnReport  	    Alignment report to be updated.  See <http://www.mothur.org/wiki/Remove.seqs#alignreport_option>
        contigsReport	Congtigs report to be updated.  See <http://www.mothur.org/wiki/Screen.seqs#contigsreport>
        summaryFile	    Sumamry file to be updated.  See <http://www.mothur.org/wiki/Screen.seqs#summary>
    :param mustFilter:  If True, at least one filter must be supplied.
    :param mustUpdate:  If True, at least one update file must be supplied.
    :return:    An option string with the selected options.  e.g "start=100, end=300"
    """
    filters = ["start", "end", "minlength", "maxlength", "maxambig", "maxn", "maxhomop"]
    updateFiles = ["groups", "names", "alnReport", "contigsReport", "summaryFile"]
    filterString = buildCSL(args,filters)
    updateString = buildCSL(args, updateFiles)
    updateString = buildCSL(args, updateFiles)
    didntFindFilter = (len(filterString) ==  0)
    didntFindUpdate = (len(updateString) == 0)
    if mustFilter and didntFindFilter:
        raise MissingMothurFilterException

    if mustUpdate and didntFindUpdate:
        raise MissingMothurFileException

    comma = ", "
    if didntFindFilter or didntFindUpdate:
        comma = ""
    # print updateString
    # print filterString
    return "%s%s%s" % (updateString, comma, filterString)


def buildCSL(item, attributes):
    """Tests if item has the attributes in options.  Any found attributes are built into a comma-delimited string in the
            format attribute_name=attribute_value.  e.g. "attr1=5, attr2=abc, attrk=def"

    :param item:        An object to probe.
    :param attributes:  A list of the desired attributes as strings.  e.g ["attr1", "attr2", "attrk"]
    :return: A comma-delimited string in the format attribute_name=attribute_value. e.g. "attr1=5, attr2=abc, attrk=def"
    """
    optionString = ""
    for attr in attributes:
        val = ""
        try:
            val = getattr(item, attr)
            if val != None:
                optionString = optionString + "%s=%s, " % (attr, val)
        except:
            pass

    # chop off the trailing comma
    return optionString [:-2]



def splitFileBySample(fastaFile, groupsFile, splitFastaDir, memLimit=0):
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


    if memLimit:
        fileSize = os.path.getsize(fastaFile)
        if fileSize > memLimit:
                raise Exception ("Target file is too big.")

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
