from converters.capitalizeSeqs import capitalizeSeqs
from parsers.parseCROPoutToGroups import parseCROPoutToGroups
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from parsers.parseUCtoGroups import parseUCtoGroups
from renamers.renameWithoutCount import removeCountsFromGroupsFile
from renamers.updateGroups import update_groups
from classes.Helpers import *

# TODO doc
def cluster_main(input_f, outdir, program, groupsfile, threads, aux_params):
    """Switchboard for clustering options.

    :param input_f:
    :param outdir:
    :param program:
    :param groupsfile:
    :param threads:
    :param aux_params:
    :return:
    """
    # Make the output directory, failing if it already exists

    makeDirOrdie(outdir)
    # Grab the fasta file(s) to cluster_main
    inputs = getInputFiles(input_f)
    debugPrintInputInfo(inputs, "clustered")
    pool = init_pool(min(len(inputs), threads))
    printVerbose("\nClustering with: %s\n" % program)
    rslt = (None,None)
    if program == "crop":
        required = (inputs, outdir, groupsfile, pool, )
        keys = ["blocksize", "clustpct", "maxmcmc", "maxsm", "rare"]
        missing_param =False
        for key in keys:
            val = aux_params.get(key, None)
            if val is None:
                print "Error: '%s' parameter is requried for crop clustering" % key
                missing_param = True
        if not missing_param:
            params = list(required)
            for key in keys:
                params.append(aux_params.get(key, None))
            print params
            rslt = cluster_crop(*params)

    elif program == "vsearch":
        vsearch_id_pct = aux_params.get('idpct', None)
        if vsearch_id_pct is None:
            print "Error: 'idpct' parameter is required for vsearch clustering"
        rslt = cluster_vsearch(inputs, outdir, groupsfile, pool, vsearch_id_pct)

    else: #"swarm"
        rslt = cluster_swarm(inputs, outdir, groupsfile, pool)

    # Grab the resulting file lists from clustering
    aux_file_list, groups_file_list = rslt

    # Move the final groups file(s) to the groups dir
    groups_dir = makeDirOrdie("%s_groups_files" % outdir)
    bulk_move_to_dir(groups_file_list, groups_dir)

    # Move aux files to the aux dir
    aux_dir = makeAuxDir(outdir)
    bulk_move_to_dir(aux_file_list, aux_dir)

    # Cleanup the pool
    cleanup_pool(pool)


def handle_groups_file_update(outdir, groupsfile, clustering_groups_files_uncount):
    """Checks if the user specified groups file exists, and updates the groupings with clustering data.
        Returns a list of the most up to date groups files.

    :param groupsfile: Optional filepath to the .groups file, or a folder of .groups files to use as a reference.
    :param clustering_groups_files_uncount: The output groups file from clustering, with trailing replication counts
        removed from sequence names.  Names in this file should match those used in the user-specified groups file
        groupsfile.
    :return: A list of filenames pointing to the most up to date groups files.
    """
    most_recent_groups_files = clustering_groups_files_uncount
    if groupsfile:
        # Try to grab groups files
        user_specified_groups_files = getInputFiles(groupsfile, critical=False)
        # If we have files at the given location
        if len(user_specified_groups_files) != 0:
            most_recent_groups_files = user_specified_groups_files
            printVerbose("Updating .groups files with clustering data")
            debugPrintInputInfo(most_recent_groups_files, "used as groups references")
            update_groups(most_recent_groups_files, clustering_groups_files_uncount, outdir, "postcluster")
            printVerbose("Done updating .groups files.")
            most_recent_groups_files = getInputFiles(outdir, "postcluster*.groups")
    else:
        printVerbose("No name files provided, assuming singletons...\n")
    return most_recent_groups_files


def cluster_vsearch(inputs, outdir, groupsfile, pool, idpct):
    """Clusters sequences using SWARM.
    :param : An argparse object with the following parameters:
                    input       A file or folder containing fasta files to cluster_main.
                    output      The output directory results will be written to.
                    groupsfile	A groups file or folder containinggroups files that describe the input.
                                Note: if no groups file is supplied, then entries in the fasta file are assumed to be
                                    singleton sequences.
                    idpct       Real number in the range (0,1] that specifies the minimum simmilarity threshold for
                                    clustering.  e.g. .95 indicates that a candidate sequence 95% must be at least
                                    95% simmilar to the seed sequence to be included in the cluster_main.
    """
    # RUN CLUSTERING
    # " --cluster_size %s -id %f --centroids %s  --uc %s",
    parallel(runProgramRunnerInstance,
             [ProgramRunner(ProgramRunnerCommands.CLUSTER_VSEARCH, [input_, float(idpct),
                                                                    "%s/%s_seeds.fasta" % (outdir,strip_ixes(input_)),
                                                                    "%s/%s_clustered_uc" % (outdir,strip_ixes(input_)),],
                                                      {"exists": [input_]}) for input_ in inputs], pool)

    # PARSE UC FILE TO GROUPS FILE
    printVerbose("Parsing the clustered uc files to groups files")
    clustered_uc_files = getInputFiles(outdir, "*_clustered_uc")
    debugPrintInputInfo(clustered_uc_files, "parsed to groups")
    parallel(runPythonInstance,
             [(parseUCtoGroups, input_, "%s/%s.groups" % (outdir, strip_ixes(input_)))
              for input_ in clustered_uc_files], pool)

    # REMOVE COUNTS FROM CLUSTERING GROUPS FILE
    printVerbose("Cleaning the .groups file from clustering")
    # Grab the current groups file and the new clustered groups file (which needs to be cleaned)
    clustered_groups_files = getInputFiles(outdir, "*_clustered.groups")
    # Remove counts from the clustering groups files
    debugPrintInputInfo(clustered_groups_files, "cleaned")
    parallel(runPythonInstance,
             [(removeCountsFromGroupsFile, input_, "%s/%s_uncount.groups" % (outdir, strip_ixes(input_)))
              for input_ in clustered_groups_files], pool)
    printVerbose("Done cleaning groups files.")

    # Collect the groups file from clustering with counts removed
    cleaned_clustered_groups_files = getInputFiles(outdir, "*_uncount.groups")

    # Resolve the user specified names file if necessary
    final_groups_files = handle_groups_file_update(outdir, groupsfile, cleaned_clustered_groups_files)

    # GATHER AUX FILES
    aux_files = getInputFiles(outdir, "*", "*_seeds.fasta")

    return aux_files, final_groups_files


def cluster_crop(inputs, outdir, groupsfile, pool, blocksize, clustpct, maxmcmc, maxsm, rare):
    """Clusters sequences using CROP.
    :param : An argparse object with the following parameters:
                    input       A file or folder containing fasta files to cluster_main.
                    output      The output directory results will be written to.
                    groupsfile	A groups file or folder containinggroups files that describe the input.
                                Note: if no groups file is supplied, then entries in the fasta file are assumed to be
                                    singleton sequences.
    """

    blockcount = ""
    if blockcount:
        blockcount = "-b %d" % blockcount

    # RUN CLUSTERING
    # crop -i %s -o %s -z %s -c %s -e %s -m %s%s
    parallel(runProgramRunnerInstance,
             [ProgramRunner(ProgramRunnerCommands.CLUSTER_CROP,
                            [input_, "%s/%s" % (outdir, strip_ixes(input_)), blocksize, clustpct,
                                maxmcmc, maxsm, rare, blockcount],
                            {"exists": [input_]}) for input_ in inputs], pool)

    # CLEAN THE OUTPUT GROUPS FILE
    printVerbose("Parsing the groups file from clustering")
    clustered_groups_files = getInputFiles(outdir, "*.cluster_main.list")
    debugPrintInputInfo(clustered_groups_files, "converted to groups files")
    parallel(runPythonInstance,
             [(parseCROPoutToGroups, input_, "%s/%s_uncount.groups" % (outdir, strip_ixes(input_)))
              for input_ in clustered_groups_files], pool)
    printVerbose("Done parsing groups file.")

    # Collect the groups file from clustering with counts removed
    cleaned_clustered_groups_files = getInputFiles(outdir, "*_uncount.groups")

    # Resolve the user specified names file if necessary
    final_groups_files = handle_groups_file_update(outdir, groupsfile, cleaned_clustered_groups_files)

    # GATHER AUX FILES
    input_dir = getDirName(input)
    aux_files = cleaned_clustered_groups_files
    aux_files += getInputFiles(input_dir, "*.unique")
    aux_files += getInputFiles(input_dir, "*.unique.list")
    aux_files += getInputFiles(input_dir, "*.unique.TempCenters.Rare")
    aux_files += getInputFiles(outdir, "*.cluster_main")
    aux_files += getInputFiles(outdir, "*.cluster_main.list")
    aux_files += getInputFiles(outdir, "*.log")
    aux_files += getInputFiles(".", "LikelihoodRatio.txt")
    return aux_files, final_groups_files


# TODO IMPORTANT: make sure that the abundances in dereplicated_renamed.fasta are sorted in decreasing order
# We need to explore this more. Run by default for now....
# Nore that any program can be used here as long as two files
# 1- seeds: contains the cluster_main centroids. This file contains the updates counts for each cluster_main.
# ex. a seq 97_2 from the cluster_main, if selected as a seed, would not be for example 97_100. This indicates
# that 98 sequences are now assigned to the cluster_main for which 97 is a seed.
# 2-clustering.out: contains the clustering results. (see file for sample format)

# ~/bin/swarm/src/swarm dereplicated_renamed.fasta \
#  		      				       --output-file clustering.out -u uclust_file -w seeds
# "swarm": program_paths["SWARM"] + " \"%s\" --output-file \"%s\" \
#                                            -u \"%s\" -w \"%s\"",

def cluster_swarm(inputs, outdir, groupsfile, pool):
    """Clusters sequences using SWARM.
    :param : An argparse object with the following parameters:
                    input       A file or folder containing fasta files to cluster_main.
                    output      The output directory results will be written to.
                    groupsfile	A groups file or folder containinggroups files that describe the input.
                                Note: if no groups file is supplied, then entries in the fasta file are assumed to be
                                    singleton sequences.
    """

    # RUN CLUSTERING
    parallel(runProgramRunnerInstance,
             [ProgramRunner(ProgramRunnerCommands.CLUSTER_SWARM,
                            [input_, "%s/%s_clustered" % (outdir, strip_ixes(input_)),
                                   "%s/%s_clustered_uc" % (outdir, strip_ixes(input_)),
                                   "%s/%s_clustered_seeds" % (outdir, strip_ixes(input_))],
                              {"exists": [input_]}) for input_ in inputs], pool)

    # PARSE UC FILE TO GROUPS FILE
    printVerbose("Parsing the clustered uc files to groups files")
    clustered_uc_files = getInputFiles(outdir, "*_clustered_uc")
    debugPrintInputInfo(clustered_uc_files, "parsed to groups")
    parallel(runPythonInstance,
             [(parseUCtoGroups, input_, "%s/%s.groups" % (outdir, strip_ixes(input_)))
              for input_ in clustered_uc_files], pool)
    printVerbose("Done parsing groups files.")

    # REMOVE COUNTS FROM CLUSTERING GROUPS FILE
    printVerbose("Cleaning the .groups file from clustering")
    # Grab the current groups file and the new clustered groups file (which needs to be cleaned)
    clustered_groups_files = getInputFiles(outdir, "*_clustered.groups")
    debugPrintInputInfo(clustered_groups_files, "cleaned")
    parallel(runPythonInstance,
             [(removeCountsFromGroupsFile, input_, "%s/%s_uncount.groups" % (outdir, strip_ixes(input_)))
              for input_ in clustered_groups_files], pool)
    printVerbose("Done cleaning groups files.")

    printVerbose("Capitalizing sequences")
    # Convert the seeds files to uppercase (swarm writes in lowercase)
    inputs = getInputFiles(outdir, "*_seeds")
    parallel(runPythonInstance,
             [(capitalizeSeqs, input_, "%s.fasta" % input_) for input_ in inputs], pool)
    printVerbose("Done capitalizing sequences")

    # Collect the groups file from clustering with counts removed
    cleaned_clustered_groups_files = getInputFiles(outdir, "*_uncount.groups")

    # Resolve the user specified names file if necessary
    final_groups_files = handle_groups_file_update(outdir, groupsfile, cleaned_clustered_groups_files)

    # GATHER AUX FILES
    aux_files = getInputFiles(outdir, "*", "*_seeds.fasta")
    return aux_files, final_groups_files

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Usage: file_to_clean    output_file    list_of_chars_to_remove   filetype "
    else:
        cluster_main(*sys.argv[1:])
