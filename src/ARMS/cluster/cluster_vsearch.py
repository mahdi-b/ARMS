from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from cluster_helpers import handle_groups_file_update
from parse.parseUCtoGroups import parseUCtoGroups
from rename.renameWithoutCount import removeCountsFromGroupsFile


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

