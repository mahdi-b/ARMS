import os
import sys
from countToDict import parseCountFileToDict
from utils.joinFiles import joinFiles

def updateNames(old_names_files, new_names_files, out_dir):
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


    :return: Filepath to the updated names file
    """
    old_names_temp_file = "%s/temp_old_merged.names" % out_dir
    new_names_temp_file = "%s/temp_new_merged.names" % out_dir
    output_name = "%s/UPDATED.names" % out_dir

    # Concat the old and new names files respectively
    joinFiles(old_names_files, old_names_temp_file)
    joinFiles(new_names_files, new_names_temp_file)

    # parse the files to dictionaries
    old_seeds = parseCountFileToDict(old_names_temp_file)
    new_seeds = parseCountFileToDict(new_names_temp_file)


    old_keys = old_seeds.keys()
    with open(output_name, 'w') as output:
        for new_seed in new_seeds.keys():
            entries_in_new_seed = new_seeds[new_seed]
            temp = ""
            for entry in entries_in_new_seed.split(' '):
                if (entry in old_keys):
                    temp += " " + old_seeds[entry]
            print temp
            output.write("%s\t%s\n" % (new_seed, entries_in_new_seed + temp))
    output.close()
    os.remove(old_names_temp_file)
    os.remove(new_names_temp_file)
    return output_name

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fasta count_file outfile"
        exit()
    updateNames(*sys.argv[1:4])