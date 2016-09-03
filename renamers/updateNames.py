import os
import sys
from parseNamesFileToDict import parseNamesFileToDictOfChilden
from utils.joinFiles import joinFiles

def updateNames(old_names_files, new_names_files, out_dir, out_prefix):
    """Updates an old_names file with the results of a new_names file, and writes the results to out_names.
    E.g. Given the old_names file lists:
    old_names:
    1   2 3
    4   5 6
    and the new_names file lists
    new_names:
    1   4
    Then return an out_names file listing
    out_names
    1   2 3 4 5 6

    Finer points:
    1. The list of child sequences following the seed should not contain the seed.
    2. The size of the cluster represented by the seed is the number of children succeeding the seed,
            plus one for the seed.
    :param new_names: The current iteration of the names file.
    :param old_names: The previous iteration of the names file.
    :param out_names: The resulting updated names file.
    :param out_prefix: The prefix for the output filename.


    :return: Filepath to the updated names file
    """
    if not (len(old_names_files)  and len(new_names_files)):
        print "Received empty file lists.  Aborting."
        exit()
    print "Using %s and %s to generate updated names file %s_updated.names" % \
          (old_names_files[0], new_names_files[0], out_prefix)
    old_names_temp_file = "%s/temp_old_merged.names" % out_dir
    new_names_temp_file = "%s/temp_new_merged.names" % out_dir
    output_name = "%s/%s_updated.names" % (out_dir, out_prefix)

    # Concat the old and new names files respectively
    joinFiles(old_names_files, old_names_temp_file)
    joinFiles(new_names_files, new_names_temp_file)

    # parse the names files to dictionaries of children
    old_seeds = parseNamesFileToDictOfChilden(old_names_temp_file)
    new_seeds = parseNamesFileToDictOfChilden(new_names_temp_file)


    old_keys = old_seeds.keys()
    with open(output_name, 'w') as output:
        for new_seed in new_seeds.keys():
            my_old_children = []
            children_of_my_new_children = []
            # Back in my day, I used to be a seed!
            if new_seed in old_keys:
                my_old_children = old_seeds[new_seed].split(" ")
            my_new_children = new_seeds[new_seed].split(" ")
            for entry in my_new_children:
                if (entry in old_keys):
                    children_of_my_new_children += old_seeds[entry].split(" ")
            all_my_children = list(set(my_old_children + my_new_children + children_of_my_new_children))
            output.write("%s\t%s\n" % (new_seed, " ".join(all_my_children)))
            # print "Seed %s has %d children" % (new_seed, len(all_my_children))
            print "%s_%d = %d old  +  %d new  +   %d children of new" % (new_seed, len(all_my_children),
                                         len(my_old_children), len(my_new_children), len(children_of_my_new_children))
    output.close()
    #os.remove(old_names_temp_file)
    #os.remove(new_names_temp_file)
    return output_name

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fasta count_file outfile"
    else:
        updateNames(*sys.argv[1:4])