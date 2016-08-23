import glob
import logging
import os
import re
import sys
import Validator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from multiprocessing import Pool


class printVerbose(object):
    """A class to toggle verbose printing."""

    # If True, verbose printing is enabled
    VERBOSE = False

    def __init__(self, msg, newline=True, out=sys.stdout):
        """Prints the msg to console if printVerbose.VERBOSE is True.
        :param msg: The thing to print.
        :param newline: Should a newline character be printed before the msg? Default True.
        :param out: Where to write the verbose message.  Default sys.stdout
        """
        if printVerbose.VERBOSE:
            if newline:
                out.write("\n")
            out.write(msg)


class MissingMothurFilterException(Exception):
    def __str__(self):
        return "Expected a mothur filter parameter, but found none"


class MissingMothurFileException(Exception):
    def __str__(self):
        return "Expected a mothur file to update, but found none."


def helpValidate(conditions):
    """Validates all conditions using a Validator object, returning True if successful and raising an exception
        if not.

    :param conditions: A dictionary mapping <Validator>.function names to a list of parameters.  The dictionary
                            key/function name is called on each parameter in the corresponding list of parameters.
                            All parameters must satisfy their <Validator>.function in order for the specified
                            program to execute.
    :return: True if validation is successful, and raising an exception in <Validator.function> if not.
    """
    try:
        for condition in conditions.iteritems():
            getattr(Validator, condition[0])(condition[1])
        # TODO sanitize inputs
        return True
    except Exception, error:
        print error

def runProgramRunnerInstance(my_instance):
    """Runs an instance of a ProgramRunner.  Calls ProgramRunner.run() for a ProgramRunner object.'
        :param my_instance A fully initalized ProgramRunner object to run.
    """
    # logging.info(myInstance.dryRun())
    return my_instance.run()


def runPythonInstance(params):
    """Wraps a multi-parameter, local, python method, so that it can be run as a multiprocessing.Pool job.

    :param params:  A tuple, where the first item is a function, and the remainder is a set of parameters
    :return:        The output of params[0](*params[1:]
    """
    printVerbose(str(params))
    func = params[0]
    args = params[1:]
    return func(*args)


def parallel(function, data, pool=Pool(processes=1), debug=False):
    """Executes jobs in parallel.
    :param function:    The function to call.  Generally runInstance() for ProgramRunner objects, or a local python
                            function.
    :param data: A list of arguments to run the function over.
    :param pool: An initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    :return: A list of results.
    """
    if debug:
        return data
    else:
        pool.map_async(function, data).get(999999999)



def makeDirOrdie(dir_path, orDie=True):
    """Creates a directory 'dirPath' or exits if the 'dirPath' directory already exists.  Prevents unnecessary execution.
    :param dir_path: The filepath to the directory to look for.

    :return: The path to the new directory (same as 'dir_path').
    """
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    else:
        if orDie:
            sys.exit("Directory %s already exists " % dir_path)
    return dir_path

def makeAuxDir(dir_path):
    """Create a directory 'dirPath_aux' if it doesn't already exist.
    :param dir_path: The filepath to the directory to look for.

    :return: The path to the new directory (same as 'dirPath').
    """
    aux_path = dir_path + "_aux"
    makeDirOrdie(aux_path)
    return aux_path


def cleanupPool(pool):
    """Cleans up a pool object.  Used when the user provides a KeyboardInterrupt

    :param pool: A fully initalized multiprocessing.Pool object.
    """
    pool.terminate()
    pool.join()
    exit()


def strip_ixes(path):
    """Strips prefix/suffixes from a file to get a base identifier for file operations.

    :param path: The path to clean.
    :return: A cleaned file identifier.
    """
    file_name = getFileName(path)
    # name = re.sub(r'_splitOut_\d+', '', file_name)
    # name = re.sub(r'_part_\d+', '', name)
    name = file_name
    ixes=[ "_renamed", "_debarcoded", ".assembled", ".discarded", ".unassembled", "_cleaned", "_derepCount","_derep",
           "_uc", "_splitOut", ".denovo.uchime", "_derepCount", "_uncount", "_counts", "_seeds"]
    for ix in ixes:
        name = name.replace(ix, "")
    return name

def clip_count(filename, delim=''):
    return re.sub(r'_\d+', '', filename)

def enumerateDir(dir_, pattern="*"):
    """Returns all the files in a directory, optionally requiring a pattern match.

    :param dir_:     The directory to look in/at.
    :param pattern: An optional pattern to match
    :return:
    """
    hits = "%s/%s" % (dir_, pattern)
    return sorted(glob.glob(hits), key=str.lower)


def getInputs(path, match_pattern="*", mismatch_pattern="", critical=True, ignore_empty_files=True):
    """Checks if a path is a file or folder.  Optionally takes a match and mismatch pattern.
        Returns a list of either a single file, or all files in a folder (matching match_pattern, if it was provided,
        and not matching mismatch_pattern, if it was provided).  If critical is True, throws an error if no files were
        found.  If ignoreEmptyFiles is True, empty files are not returned.

    :param path:             File path to a folder or file to look in/at.
    :param match_pattern:            Optional pattern string. e.g. *.txt grabs all .txt files in path.
    :param mismatch_pattern:         Mismatch pattern.  Ignores all files matching this pattern, even if they
                                previously matched.
    :param critical:         If True, throws an exception if no files are found.  Default True.
    :param ignore_empty_files: If True, does not return empty files.  Default True.
    :return:            A list of one or more files.
    """
    rslt = []
    if os.path.isfile(path):
        rslt = [path]
    elif os.path.isdir(path):
        patternMatch = enumerateDir(path, match_pattern)
        if mismatch_pattern:
            negative_match = enumerateDir(path, mismatch_pattern)
            rslt = list(set(patternMatch) - set(negative_match))
        else:
            rslt = patternMatch

    else:
        logging.error("Found no  matching inputs matching %s at %s" % (match_pattern, path))
        if critical:
            print "Error: Found no matching inputs matching %s at %s" % (match_pattern, path)
            exit()
    if ignore_empty_files:
        rslt = [path for path in rslt if os.path.getsize(path)]
    return rslt


def getDir(path):
    """Returns the last directory name in path.

    :param path: The filepath.
    :return:    The last directory on path, or path if it is a directory.
    """
    if os.path.isfile(path):
        return os.path.dirname(path)
    elif os.path.isdir(path):
        return path
    else:
        print "Error: Invalid path."
        raise Exception


def getFileName(path):
    """Returns the filename (no extension, no directory) in an absoloute filepath.

    :param path:The file path.
    :return:    Returns the filename on a path with no directory prefix or file extension.
    """
    return os.path.splitext(os.path.basename(path))[0]


#TODO work on this
def sanitize_string(input_):
    """Returns a white-list sanitized string.

    :param input_:
    :return:
    """
    return input_


def sanitize(inputs):
    """White list sanitizes a dictionary of inputs, returning a dictionary of sanitized strings.  Raises an exception if
        an invalid character is found.
    :param inputs: The list of inputs to sanitize

    :return: A list of sanitized inputs
    """

    for name, input_ in inputs.iteritems():
        text = str(input_)

    return inputs


def bulk_move_to_dir(target_files, dest_dir_path):
    """Moves a list of files with a given filepath to a given directory.

    :param target_files:    A list of files to move
    :param dest_dir_path:   The directory to move the files to.
    :return:
    """
    if target_files is not None:
        for target_file in target_files:
            move_to_dir(target_file, dest_dir_path)


def move_to_dir(target_file_path, dest_dir_path):
    """Moves a file with a given filepath to a given directory.  Sanitizes inputs before execution.

    :param target_file_path:          File path to the file to be moved.  e.g. "dirA/a.txt"
    :param dest_dir_path:     File path to the file's new destination.  e.g. "dirB/"
    """
    target_file = target_file_path.split('/')[-1]
    dest_file_path = "%s/%s" % (dest_dir_path, target_file)
    move(target_file_path, dest_file_path)


def move(origin, destination):
    """Moves a file from an origin to a destination.  Sanitizes inputs before execution.

    :param origin: File path to the file to be moved.  e.g. "dirA/a.txt"
    :param destination: File path to the file's new destination.  e.g. "dirB/a.txt"
    """
    os.rename(sanitize_string(origin), sanitize_string(destination))


def debugPrintInputInfo(input_, action_suffix):
    """Prints debug information about an input.

    :param input_: The input to log debug info for.
    :param action_suffix: The past-tense operation to be preformed on the input.
    :return:
    """
    logging.debug("%d files to be %s:" % (len(input_), action_suffix))
    logging.debug(str(input_))


def mothur_parseInputFileType(args):
    """Gets the input file format and input file from the args argparser.

    :param args:
    :return: A tuple: (input file format, input file path)
    """
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
    filters = ["start", "end", "minlength", "maxlength", "maxambig", "maxn", "maxhomop", "pdiff", "bdiff"]
    update_files = ["groups", "names", "alnReport", "contigsReport", "summaryFile"]
    filter_string = buildCSL(args, filters)
    update_string = buildCSL(args, update_files)
    didnt_find_filter = (len(filter_string) == 0)
    didnt_find_update = (len(update_string) == 0)
    if mustFilter and didnt_find_filter:
        raise MissingMothurFilterException

    if mustUpdate and didnt_find_update:
        raise MissingMothurFileException

    comma = ", "
    if didnt_find_filter or didnt_find_update:
        comma = ""
    # print update_string
    # print filter_string
    return "%s%s%s" % (update_string, comma, filter_string)


def buildCSL(item, attributes):
    """Tests if item has the attributes in options.  Any found attributes are built into a comma-delimited string in the
            format attribute_name=attribute_value.  e.g. "attr1=5, attr2=abc, attrk=def"

    :param item: An object to probe.
    :param attributes:  A list of the desired attributes as strings.  e.g ["attr1", "attr2", "attrk"]
    :return: A comma-delimited string in the format attribute_name=attribute_value. e.g. "attr1=5, attr2=abc, attrk=def"
    """
    option_string = ""
    for attr in attributes:
        try:
            val = getattr(item, attr)
            if val is not None:
                option_string = "%s%s=%s, " % (option_string, attr, val)
        except:
            pass

    # chop off the trailing comma
    return option_string[:-2]



def splitFileBySample(fastaFile, groupsFile, splitFastaDir, mem_limit=0):
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

    if mem_limit:
        file_size = os.path.getsize(fastaFile)
        if file_size > mem_limit:
                raise Exception("Target file is too big.")

    # TODO: test that the file fits into memory, otherwise this could cause problems
    my_sequences = SeqIO.to_dict(SeqIO.parse(fastaFile, 'fasta'))

    seqs_by_group = defaultdict(list)
    for line in open(groupsFile, 'r'):
        seq, group = line.rstrip().split()
        seqs_by_group[group].append(seq)
    for group in seqs_by_group.keys():
        out_file = open(os.path.join(splitFastaDir, group), 'a')
        for seq in seqs_by_group[group]:
            SeqIO.write(
                SeqRecord(my_sequences[seq].seq.ungap(".").ungap("-"), id=seq, description=""), out_file, 'fasta')
        out_file.close()
