from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *
from parse.parseUCtoGroups import parseUCtoGroups
from rename.renameWithoutCount import removeCountsFromGroupsFile
from util.capitalizeSeqs import capitalizeSeqs

from classes.Helpers import *
from cluster_helpers import handle_groups_file_update


class Cluster_Program_Swarm(ChewbaccaProgram):
    name = "swarm"


    def execute_program(self):
        args = self.args
        self.cluster_swarm(args.input_f, args.outdir, args.groupsfile, args.threads)


    # TODO IMPORTANT: make sure that the abundances in dereplicated_renamed.fasta are sorted in decreasing order
    # We need to explore this more. Run by default for now....
    # Nore that any program can be used here as long as two files
    # 1- seeds: contains the cluster centroids. This file contains the updates counts for each cluster.
    # ex. a seq 97_2 from the cluster, if selected as a seed, would not be for example 97_100. This indicates
    # that 98 sequences are now assigned to the cluster for which 97 is a seed.
    # 2-clustering.out: contains the clustering results. (see file for sample format)

    # ~/bin/swarm/src/swarm dereplicated_renamed.fasta \

    def cluster_swarm(self, input_f, outdir, groupsfile, threads):
        """Clusters sequences using SWARM.
        :param : An argparse object with the following parameters:
                        input       A file or folder containing fasta files to cluster.
                        output      The output directory results will be written to.
                        groupsfile	A groups file or folder containinggroups files that describe the input.
                                    Note: if no groups file is supplied, then entries in the fasta file are assumed to be
                                        singleton sequences.
        """
        makeDirOrdie(outdir)
        # Grab the fasta file(s) to cluster
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "clustered")
        pool = init_pool(min(len(inputs), threads))

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

        # Move the final groups file(s) to the groups dir
        groups_dir = makeDirOrdie("%s_groups_files" % outdir)
        bulk_move_to_dir(final_groups_files, groups_dir)

        # Move aux files to the aux dir
        aux_dir = makeAuxDir(outdir)
        bulk_move_to_dir(aux_files, aux_dir)

        # Cleanup the pool
        cleanup_pool(pool)