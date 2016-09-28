from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from cluster_helpers import handle_groups_file_update
from parse.parseCROPoutToGroups import parseCROPoutToGroups

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