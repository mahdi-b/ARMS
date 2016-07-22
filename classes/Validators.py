import os


def exists(files):
    """Returns True if all files in a list exist in the current directory, and raises an Exception if not.

    :param files: A list of files to test
    :return: True if all files exist, raises an Exception if not.
    """
    for file_ in files:
        if not os.path.exists(file):
            raise Exception("File %s does not exist, executions halted" % file_)
    return True


def installed(progName):
    """Tests if a program is installed, and available on $PATH.

    :param progName: Command line monicker of the program to test.
    :return: True if all programs exist, and raises an Exception if not.
    """
    pass
