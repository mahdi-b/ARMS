
# parses a count file and returns a dict mapping each unique seed to its count
def parseCountFileToCountDict(count_file):
    seeds = {}
    print "Reading count file: %s" % count_file
    # collect the seed names, and the children sequence names
    i = 0
    nb_lines = 0
    for line in open(count_file, 'r'):
        seed, children = line.rstrip().split("\t")
        print line
        seeds[seed] = len(children.split(" "))
        if nb_lines % 1000000 == 0:
            print "%s lines processed" % nb_lines
        nb_lines +=1
    print "Done reading count file."
    return seeds

def parseCountFileToDict(count_file):
    seeds = {}
    print "Reading count file: %s" % count_file
    # collect the seed names, and the children sequence names
    i = 0
    nb_lines = 0
    for line in open(count_file, 'r'):
        seed, children = line.rstrip().split("\t")
        print line
        seeds[seed] = children
        if nb_lines % 1000000 == 0:
            print "%s lines processed" % nb_lines
        nb_lines +=1
    print "Done reading count file."
    return seeds