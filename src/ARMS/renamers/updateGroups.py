import os
import sys
from utils.joinFiles import joinFiles
from parsers.parseGroupsFileToDict import parseGroupsFileToDictOfChilden


def update_groups(old_groups_files, new_groups_files, out_dir, out_prefix):
    """Updates an old_groups file with the results of a new_groups file, and writes the results to a new groups file.
    E.g. Given the old_groups file lists:
    old_groups:
    1   2 3
    4   5 6
    and the new_groups file lists
    new_groups:
    1   4
    Then return an out_groups file listing
    out_groups
    1   2 3 4 5 6

    Finer points:
    1. The list of child sequences following the seed should not contain the seed.
    2. The size of the cluster_main represented by the seed is the number of children succeeding the seed,
            plus one for the seed.
    :param new_groups_files: The current iteration of the groups file.
    :param old_groups_files: The previous iteration of the groups file.
    :param out_dir: The resulting updated groups file.
    :param out_prefix: The prefix for the output filename.


    :return: Filepath to the updated groups file
    """
    if not (len(old_groups_files) and len(new_groups_files)):
        print "Received empty file lists.  Aborting."
        return
    print "Using %s and %s to generate updated groups file %s_updated.groups" % \
          (old_groups_files[0], new_groups_files[0], out_prefix)
    old_groups_temp_file = "%s/temp_old_merged.groups" % out_dir
    new_groups_temp_file = "%s/temp_new_merged.groups" % out_dir
    output_file = "%s/%s_updated.groups" % (out_dir, out_prefix)

    # Concat the old and new groups files respectively
    joinFiles(old_groups_files, old_groups_temp_file)
    joinFiles(new_groups_files, new_groups_temp_file)

    # parse the groups files to dictionaries of children
    old_seeds = parseGroupsFileToDictOfChilden(old_groups_temp_file)
    new_seeds = parseGroupsFileToDictOfChilden(new_groups_temp_file)
    new_keys = new_seeds.keys()
    total = len(new_keys)
    i = 0
    outstring = ""
    with open(output_file, 'w') as output:
        for new_seed in new_keys:
            i += 1
            if i % 100000 == 0:
                print "Processed %d / %d lines\n" % (i, total)
                output.write(outstring)
                outstring = ""
            my_old_children = []
            children_of_my_new_children = []
            # Back in my day, I used to be a seed!
            if old_seeds.has_key(new_seed):
                my_old_children = old_seeds[new_seed].split(" ")
            my_new_children = new_seeds[new_seed].split(" ")
            for entry in my_new_children:
                if old_seeds.has_key(entry):
                    children_of_my_new_children += old_seeds[entry].split(" ")
            # list(set(my_old_children + my_new_children + children_of_my_new_children))
            all_my_children = my_old_children + my_new_children + children_of_my_new_children

            outstring += "%s\t%s\n" % (new_seed, " ".join(set(all_my_children)))
            # print "Seed %s has %d children" % (new_seed, len(all_my_children))
            # print("%s_%d = %d old  +  %d new  +   %d children of new" %
            #           (new_seed, len(all_my_children), len(my_old_children), len(my_new_children),
            #            len(children_of_my_new_children)))
        output.write(outstring)
        output.close()
    os.remove(old_groups_temp_file)
    os.remove(new_groups_temp_file)
    return output_file


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Usage: old_groups_files, new_groups_files, out_dir, out_prefix"
    else:
        update_groups(*sys.argv[1:5])
