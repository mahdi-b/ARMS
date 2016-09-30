from classes.ChewbaccaProgram import *
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from parse.parseCROPoutToGroups import parseCROPoutToGroups

from classes.Helpers import *
from cluster_helpers import handle_groups_file_update


class Cluster_Program_Crop(ChewbaccaProgram):
    """Uses CROP to cluster a set of input sequences.
    """
    name = "crop"


    def execute_program(self):
        args = self.args

        blockcount = ""
        if args.blockcount:
            blockcount = "-b %d" % args.blockcount

        self.cluster_crop(args.input_f, args.outdir, args.groupsfile, args.threads,
                          args.blocksize, args.clustpct, args.maxmcmc, args.maxsm, args.rare, blockcount)


    def cluster_crop(self, input_f, outdir, groupsfile, threads, blocksize, clustpct, maxmcmc, maxsm, rare, blockcount):
        """Clusters sequences using CROP.

        :param input_f: Filepath to the input fasta file to cluster.
        :param outdir: Filepath to the output directory.
        :param groupsfile: Filepath to the groups file to use as a reference for dereplication counting.
        :param threads: The maximum number of processes to use to cluster.
        :param blocksize: Size of blocks to be used for all rounds (if -b is specified, then -z will not affect the
                            first round.  For data set with different average sequence length, this parameter should \
                            be tuned such that it won't take too long for each block to do pariwise alignment.  Hint \
                            for choosing z: z*L<150,000, where L is the average length of the sequences.
        :param clustpct: The minimum similarity threshold for clustering.  Either 'g' for 95% or 's' for 97%.
        :param maxmcmc: This parameter specifies the number of iterations of MCMC. Default value is 2000. Increase \
                            this value to enhance accuracy (recommended value is at least 10*block size).
        :param maxsm: This parameter specifies the maximum number of 'split and merge' process to run.  Max is 20.
        :param rare: The maximum cluster size allowed to be classified as 'rare'. Clusters are defined as either \
                            'abundant' or 'rare'. 'Abundant' clusters will be clustered first, then the 'rare' \
                            clusters are mapped to the 'abundant' clusters.  Finally, 'rare' clusters which cannot be \
                            mapped will be clustered separately. e.g. If r=5, the clusters with size <=5 will be \
                            considered 'rare' in above procedure. and r=0 will yield the best accuracy. If you \
                            believe your data is not too diverse to be handled, then r=0 will be the best choice.
        """

        makeDirOrdie(outdir)
        # Grab the fasta file(s) to cluster
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "clustered")
        pool = init_pool(min(len(inputs), threads))



        # RUN CLUSTERING
        # crop -i %s -o %s -z %s -c %s -e %s -m %s%s
        parallel(runProgramRunnerInstance,
                 [ProgramRunner(ProgramRunnerCommands.CLUSTER_CROP,
                                [input_, "%s/%s" % (outdir, strip_ixes(input_)), blocksize, clustpct,
                                    maxmcmc, maxsm, rare, blockcount],
                                {"exists": [input_]}) for input_ in inputs], pool)

        # CLEAN THE OUTPUT GROUPS FILE
        printVerbose("Parsing the groups file from clustering")
        clustered_groups_files = getInputFiles(outdir, "*.cluster.list")
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
        input_dir = getDirName(input_f)
        aux_files = cleaned_clustered_groups_files
        aux_files += getInputFiles(input_dir, "*.unique")
        aux_files += getInputFiles(input_dir, "*.unique.list")
        aux_files += getInputFiles(input_dir, "*.unique.TempCenters.Rare")
        aux_files += getInputFiles(outdir, "*.cluster")
        aux_files += getInputFiles(outdir, "*.cluster.list")
        aux_files += getInputFiles(outdir, "*.log")
        aux_files += getInputFiles(".", "LikelihoodRatio.txt")

        # Move the final groups file(s) to the groups dir
        groups_dir = makeDirOrdie("%s_groups_files" % outdir)
        bulk_move_to_dir(final_groups_files, groups_dir)

        # Move aux files to the aux dir
        aux_dir = makeAuxDir(outdir)
        bulk_move_to_dir(aux_files, aux_dir)

        # Cleanup the pool
        cleanup_pool(pool)