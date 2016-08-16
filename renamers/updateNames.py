import sys
from countToDict import parseCountFileToDict

def updateNames(old_names, new_names, out_names):
    """Updates an old_names file with the results of a new_names file, and writes the results to out_names.
    E.g. Given the old_ames file lists:
    old_names:
    1   2 3
    4   5 6
    and the new_names file lists
    new_names:
    1   4
    Then return an out_names file listing
    out_names
    1   2 3 4 5 6

    :param new_names: The current iteration of the names file.
    :param old_names: The previous iteration of the names file.
    :param out_names: The resulting updated names file.
    """

    new_seeds = parseCountFileToDict(new_names)
    old_seeds = parseCountFileToDict(old_names)

    old_keys = old_seeds.keys()
    print old_seeds
    print old_keys
    with open(out_names, 'w') as output:
        for new_seed in new_seeds.keys():
            entries_in_new_seed = new_seeds[new_seed]
            temp = ""
            for entry in entries_in_new_seed.split(' '):
                if (entry in old_keys):
                    temp += " " + old_seeds[entry]
            print temp
            output.write("%s\t%s\n" % (new_seed, entries_in_new_seed + temp))
    output.close()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fasta count_file outfile"
        exit()
    updateNames(*sys.argv[1:4])