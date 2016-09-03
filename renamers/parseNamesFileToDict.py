

def parseNamesFileToDictOfCounts(names_file):
    """Given a .names file, returns a dictionary mapping each seed to the number of children it represents.

    :param names_file: A names file.
    :return: A dictionary where each seed name is a key to a count of its children.
    """
    return parseNamesFileToDict(names_file, "counts")


def parseNamesFileToDictOfChilden(names_file):
    """Given a .names file, returns a dictionary mapping each seed to a space-delimited string of its children's names.

       :param names_file: A names file.
       :return: A dictionary where each seed name is a key to a space-delimited string of its children's names.
       """
    return parseNamesFileToDict(names_file, "children")


def parseNamesFileToDict(names_file, thing_to_map):
    """Given a .names file, returns a dictionary mapping each seed to either a count of its children, or a
        space-delimited string of its children's names.

    :param names_file: A .names file.
    :param thing_to_map: Specify 'children' to map seed names to a space-delimited string of children names, or 'counts'
                            to map seed names to a count of children.
    :return: A dictionary mapping each seed to either a count of its children, or a space-delimited string of its
            children's names
    """
    seeds = {}
    print "Reading count file: %s" % names_file
    # collect the seed names, and the children sequence names
    i = 0
    nb_lines = 0
    for line in open(names_file, 'r'):
        data = line.rstrip().split("\t")
        seed = data[0]
        children = ""
        if thing_to_map == "children":
            if len(data) > 1:
                children = ' '.join(list(set(data[1:])))
            seeds[seed] = children
        if thing_to_map == "counts":
            if len(data) > 1:
                children = data[1]
            seeds[seed] = len(children.split(" "))

        if nb_lines % 1000000 == 0:
            print "%s lines processed" % nb_lines
        nb_lines +=1
    print "Done reading count file."
    return seeds