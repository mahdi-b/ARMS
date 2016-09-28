from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from cluster_helpers import handle_groups_file_update
from util.capitalizeSeqs import capitalizeSeqs
from parse.parseUCtoGroups import parseUCtoGroups
from rename.renameWithoutCount import removeCountsFromGroupsFile





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