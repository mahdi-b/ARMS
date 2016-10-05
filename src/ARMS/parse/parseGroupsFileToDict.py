from classes.Helpers import *


def parseGroupsFileToDictOfCounts(groups_file):
    """Given a .groups file, returns a dictionary mapping each seed to the number of children it represents.

    :param groups_file: A groups file.
    :return: A dictionary where each seed name is a key to a count of its children.
    """
    return parseGroupsFileToDict(groups_file, "counts")


def parseGroupsFileToDictOfChilden(groups_file):
    """Given a .groups file, returns a dictionary mapping each seed to a space-delimited string of its children's names.

       :param groups_file: A groups file.
       :return: A dictionary where each seed name is a key to a space-delimited string of its children's names.
       """
    return parseGroupsFileToDict(groups_file, "children")


def parseGroupsFileToDict(groups_file, thing_to_map):
    """Given a .groups file, returns a dictionary mapping each seed to either a count of its children, or a
        space-delimited string of its children's names.

    :param groups_file: A .groups file.
    :param thing_to_map: Specify 'children' to map seed names to a space-delimited string of children names, or 'counts'
                            to map seed names to a count of children.
    :return: A dictionary mapping each seed to either a count of its children, or a space-delimited string of its
            children's names
    """
    groups = {}
    printVerbose("Reading count file: %s" % groups_file)
    # collect the seed names, and the children sequence names
    i = 0
    nb_lines = 0
    for line in open(groups_file, 'r'):
        nb_lines +=1
        data = line.rstrip().split("\t")
        seed = data[0]
        children = ""
        if thing_to_map == "children":
            if len(data) > 1:
                children = ' '.join(list(set(data[1:])))
            groups[seed] = children
        if thing_to_map == "counts":
            if len(data) > 1:
                children = data[1]
            groups[seed] = len(children.split(" "))

        if nb_lines % 100000 == 0:
            printVerbose("%s lines processed" % nb_lines)
    printVerbose("Done reading count file.")
    return groups