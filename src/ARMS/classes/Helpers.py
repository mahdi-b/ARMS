import glob
import logging
import os
import sys
from multiprocessing import Pool
import classes.Validator
from shutil import copy2

class printVerbose(object):
    """A class to toggle verbose printing."""

    # If True, verbose printing is enabled
    VERBOSE = False

    def __init__(self, msg,out=sys.stdout):
        """Prints the msg to console if printVerbose.VERBOSE is True.
        :param msg: The thing to print.
        :param out: Where to write the verbose message.  Default sys.stdout
        """
        if printVerbose.VERBOSE:
            out.write("\n%s\n" % msg)

def helpValidate(conditions):
    """Validates all conditions using a Validator object, returning True if successful and raising an exception
        if not.

    :param conditions: A dictionary mapping <Validator>.function names to a list of parameters.  The dictionary
                            key/function name is called on each parameter in the corresponding list of parameters.
                            All parameters must satisfy their <Validator>.function in order for the specified
                            program to execute.
    :return: True if validation is successful, and raising an exception in <Validator.function> if not.
    """

    for condition in conditions.iteritems():
        getattr(classes.Validator, condition[0])(condition[1])
    return True


def runProgramRunnerInstance(my_instance):
    """Runs an instance of a ProgramRunner.  Calls ProgramRunner.run() for a ProgramRunner object.'
        :param my_instance A fully initalized ProgramRunner object to run.
    """
    # logging.info(myInstance.dry_run())
    return my_instance.run()


def runPythonInstance(params):
    """Wraps a multi-parameter, local, python method, so that it can be run as a multiprocessing.Pool job.

    :param params:  A tuple, where the first item is a function, and the remainder is a set of parameters
    :return:        The output of params[0](*params[1:]
    """
    printVerbose(str(params))
    func = params[0]
    args = params[1:]
    printVerbose("Running: %s" % str(params))
    return func(*args)


def parallel(function, data, pool):
    """Executes jobs in parallel.
    :param function:    The function to call.  Generally runInstance() for ProgramRunner objects, or a local python
                            function.
    :param data: A list of arguments to run the function over.
    :param pool: An initalized multiprocessing.Pool object.
    :return: A list of results.
    """
    try:
        return pool.map_async(function, data).get(999999999)

    except KeyboardInterrupt:
        kill_pool_and_die(pool)

    except  Exception as e:
        print e
        print sys.exc_info()
        kill_pool_and_die(pool)


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


def copy_file(target_path, dest_path):
    """Copies a file using shutil.copy2.

    :param target_path: Filepath to the file to copy.
    :param dest_path: Filepath to the destination.
    """
    copy2(target_path, dest_path)


def init_pool(requested_pool_size):
    """Initalizes a pool of a given size.

    :param requested_pool_size: The desired size of the pool
    :return: An initalized pool.
    """
    # Number of threads should be between 1 and the requested_pool_size
    # pool_size = max(1, min(requested_pool_size, cpu_count()))
    pool_size = max(1, requested_pool_size)
    printVerbose("Running with %s process(es)" % pool_size)
    pool = Pool(pool_size)
    return pool


def kill_pool_and_die(pool):
    """Cleans up a pool object.  Used when the user provides a KeyboardInterrupt

    :param pool: A fully initalized multiprocessing.Pool object.
    """
    debugPrint("\nEncountered an error!\nClosing Worker Pool...\n")
    # prevent any new jobs
    pool.close()
    # terminate the pool
    pool.terminate()
    # Wait for everything to return and exit
    pool.join()
    # New line for cursor
    print "\n"
    sys.exit()


def cleanup_pool(pool):
    """Cleans up a pool object after normal use, and exits.

    :param pool: A fully initalized multiprocessing.Pool object.
    """
    pool.close()
    pool.join()
    print "\n"
    sys.exit()


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


def clip_count(sequence_name, delim='_'):
    """Removes a trailing abundance count from a sequence name.  Used by rename sequences
        e.g. a sequence with an abundance/replication count of 10 (and a delimiter of '_'):
            "Specifically_named_seq_27_10"
            clip_count("Specifically_named_seq_27_10) -> Specifically_named_seq_27

    :param sequence_name: The sequence name with trailing replication count to process.
    :param delim: The delimiter for the trailing abundance/replication count.
    :return: The sequence name without the abundance/replication count.
    """

    return delim.join(sequence_name.split(delim)[:-1])


def enumerateDir(dir_, pattern="*"):
    """Returns all the files in a directory, optionally requiring a pattern match.

    :param dir_: The directory to look in for files.
    :param pattern: An optional pattern to match.
    :return: A of all the files in a directory.  If pattern is supplied, then return only those files matching
            the pattern.  If no files match the requested pattern, an empty list is returned.
            Returned Files are sorted by lowercase filename.
    """
    if pattern == "":
        return []
    hits = "%s/%s" % (dir_, pattern)
    return sorted(glob.glob(hits), key=str.lower)


def getInputFiles(path, match_pattern="*", mismatch_pattern="", critical=True, ignore_empty_files=True):
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


def getDirName(path):
    """If a file path is provided, returns the last dir in that path.  If a directory path is provided, returns that
        directory name."""
    if os.path.isdir(path):
        return path
    return os.path.dirname(path)


def getFileName(path):
    """Returns the filename (no trailing extension, no preceeding directory) in an absoloute filepath.

    :param path: The file path.
    :return: Returns the filename on a path with no directory prefix or file extension.
    """
    return os.path.splitext(os.path.basename(path))[0]


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
    try:
        move(target_file_path, dest_file_path)
    except:
        print "Warning: Couldn't move %s to %s\n" % (target_file, dest_dir_path)
        pass


def move(origin, destination):
    """Moves a file from an origin to a destination.  Sanitizes inputs before execution.

    :param origin: File path to the file to be moved.  e.g. "dirA/a.txt"
    :param destination: File path to the file's new destination.  e.g. "dirB/a.txt"
    """
    os.rename(origin, destination)

def debugPrint(message):
    """Logs debug messages.

    :param message: The debug message to log
    """
    logging.debug(message)


def debugPrintInputInfo(input_, action_suffix):
    """Prints debug information about an input.

    :param input_: The input to log debug info for.
    :param action_suffix: The past-tense operation to be preformed on the input.
    :return:
    """
    file_str = ""
    for f in input_:
        file_str += "\t%s\n" % str(f)

    logging.debug("%d files to be %s:\n%s\n" % (len(input_), action_suffix, file_str))


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


def validate_paired_fastq_reads(input_f, input_r):
    forwards_reads = getInputFiles(input_f, "*_forward*", critical=False)
    reverse_reads = getInputFiles(input_r, "*_reverse*", critical=False)

    if len(forwards_reads) == 0 and len(reverse_reads) == 0:
        forwards_reads = getInputFiles(input_f, "*_R1*", critical=False)
        reverse_reads = getInputFiles(input_r, "*_R2*", critical=False)

    # Ensure that we have matching left and right reads
    if len(forwards_reads) != len(reverse_reads):
        print "Error: Unequal number of forwards/reverse reads."
        return []

    if len(forwards_reads) == 0:
        print "Forwards reads should include the filename suffix \"_forward\" or \"R1\".  Reverse reads should \
                        include the filename suffix \"_reverse\" or \"R2\"."
        return
    return zip(set(forwards_reads), set(reverse_reads))


def one_to_one_or_one_to_many(n,m):
    k = len(n)
    if k ==1 or k ==len(m):
        return True
    return False


def fittosize(n,m):
    if len(n)<len(m):
        return n*len(m), m
