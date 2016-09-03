import sys
from classes.Helpers import *
# NOTE: A SEQUENCE MUST NOT APPEAR IN TWO NAMES FILES.
def buildOTUtable(latest_names_files, inital_groups_files, barcodes_file, out_file, padding=15):
    """Given a single barcodes file with all possible \
    sample names, a list of the latest names file(s), and a list of initial groups files \
    (mapping each original, undereplicated sequence to its sample name), builds an OTU \
    table and writes it to out_file.

    :param latest_names_files:  A list of the latest names files.  Any sequence name must not occur in more than one \
                                    names file.
    :param inital_groups_files: A list of the inital groups files.  This should map each sequence to its parent group.
    :param barcodes_file:       A single barcodes file listing all valid sample names.
    :param out_file:            Filepath to the output directory
    """
    seq_to_sample = {}
    max_row_name_length = 30 #TODO compute max name length
    # read the initaial groups/samples file (from rename)
    # make a single dict from all the groups/samples files mapping seqname to group
    for groups_file in inital_groups_files:
        with open(groups_file, 'r') as current_groups_file:
            for line in current_groups_file:
                name, sample = line.split("\t")
                seq_to_sample[name] = sample.rstrip()
    printVerbose("Finding sample names...")
    # read the barcodes file to get a list of all the possible groups.  sort it
    all_sample_names = []
    with open(barcodes_file, 'r') as barcodes:
        for line in barcodes:
            all_sample_names.append(line.split()[0].rstrip())
    all_sample_names.sort()
    printVerbose("Done finding keynames.  Found:")
    printVerbose(str(all_sample_names))

    with open(out_file, 'w') as out:
        header_line = "%-*s" % (max_row_name_length, "OTU")
        for sample in all_sample_names:
            header_line += "%*s" % (padding, sample)
        out.write(header_line + "\n")
        # GENERATE A DICTIONARY MAPPING SEQUENCE NAMES TO THE SAMPLE THEY CAME FROM
        # for each line in the latest names files,
        for names_file in latest_names_files:
            with open(names_file, 'r') as current_names_file:
                otu = ""
                children = ""
                # read the latest names file
                for line in current_names_file:
                    data = line.split("\t")
                    # found a cluster
                    if len(data) == 2:
                        otu = data[0].rstrip()
                        children = data[1].rstrip()

                    # found a singleton
                    elif len(data) == 1:
                        otu = data[0].rstrip()

                    # found a blank line
                    else:
                        pass

                    # GENERATE OTU ABUNDANCE BY SAMPLE
                    # initalize a count_dict with each sample as a key and a value of 0
                    sample_counts = {}
                    for sample_name in all_sample_names:
                        sample_counts[sample_name] = 0

                    # for each item in the child list:
                    for child in children.split():
                        # my_sample = lookup that item in the dict to get its sample name
                        my_sample = seq_to_sample[child]
                        # increment the abundance in that sample
                        sample_counts[my_sample] += 1
                    # don't forget the representative sequence
                    #   (the representative of the cluster is not in the child list)
                    my_sample = seq_to_sample[otu]
                    sample_counts[my_sample] += 1

                    # WRITE THE COUNTS TO THE OUT FILE
                    # for each sample in the barcodes list, write otu to a txt file as a single line
                    out_line = "%-*s" % (max_row_name_length, otu)
                    for sample_name in all_sample_names:
                        out_line += "%*s" % (padding, sample_counts[sample_name])
                    out.write(out_line + "\n")
    out.close()
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Usage: list_of_latest_names_file  list_of_inital_groups_files  barcodes_file  out_file"
        exit()
    buildOTUtable(*sys.argv[1:5])