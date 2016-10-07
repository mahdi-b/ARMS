from classes.ChewbaccaProgram import *
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from classes.Helpers import *
from Cluster_Helpers import handle_groups_file_update


class Cluster_Program_Crop(ChewbaccaProgram):
    """Uses CROP to cluster a set of input sequences.
    """
    name = "crop"

    def execute_program(self):
        args = self.args

        blockcount = ""
        if args.blockcount:
            blockcount = "-b %d" % args.blockcount

        self.cluster_crop(args.input_f, args.outdir, args.groupsfile, args.processes,
                          args.blocksize, args.clustpct, args.maxmcmc, args.maxsm, args.rare, blockcount,
                          args.extraargstring)

    def cluster_crop(self, input_f, outdir, groupsfile, processes, blocksize, clustpct, maxmcmc, maxsm, rare,
                     blockcount, extraargstring):
        """Clusters sequences using CROP.

        :param input_f: Filepath to the input fasta file to cluster.
        :param outdir: Filepath to the output directory.
        :param groupsfile: Filepath to the groups file to use as a reference for dereplication counting.
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
        :param blockcount: The size of blocks in the first round of clustering. Hint of choosing -b: Each block in the \
                            first round should contain about 50 sequences.  i.e. b=N/50, where N is the number of \
                            input sequences.  Default: # input sequences / z.
        :param processes: The maximum number of processes to use.
        :param extraargstring: Advanced program parameter string.
        """

        # Grab the fasta file(s) to cluster
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "clustered")
        pool = init_pool(min(len(inputs), processes))

        # RUN CLUSTERING
        # crop -i %s -o %s -z %s -c %s -e %s -m %s%s
        parallel(runProgramRunnerInstance,
                 [ProgramRunner(ProgramRunnerCommands.CLUSTER_CROP,
                                [input_, "%s/%s" % (outdir, strip_ixes(input_)), blocksize, clustpct,
                                    maxmcmc, maxsm, rare, blockcount],
                                {"exists": [input_]}, extraargstring) for input_ in inputs], pool)

        # CLEAN THE OUTPUT GROUPS FILE
        printVerbose("Parsing the groups file from clustering")
        clustered_groups_files = getInputFiles(outdir, "*.cluster.list")
        debugPrintInputInfo(clustered_groups_files, "converted to groups files")
        parallel(runPythonInstance,
                 [(parseCROPoutToGroups, input_, "%s/%s_uncount.groups" % (outdir, strip_ixes(input_)))
                  for input_ in clustered_groups_files], pool)
        printVerbose("Done parsing groups file.")

        # Collect the groups file from clustering with counts removed
        cleaned_clustered_groups_files = getInputFiles(outdir, "*_uncount.groups", ignore_empty_files=False)

        # Resolve the user specified names file if necessary
        final_groups_files = handle_groups_file_update(outdir, groupsfile, cleaned_clustered_groups_files)

        # GATHER AUX FILES
        input_dir = getDirName(input_f)
        aux_files = cleaned_clustered_groups_files
        aux_files += getInputFiles(input_dir, "*.unique", ignore_empty_files=False)
        aux_files += getInputFiles(input_dir, "*.unique.list", ignore_empty_files=False)
        aux_files += getInputFiles(input_dir, "*.unique.TempCenters.Rare", ignore_empty_files=False)
        aux_files += getInputFiles(outdir, "*.cluster", ignore_empty_files=False)
        aux_files += getInputFiles(outdir, "*.cluster.list", ignore_empty_files=False)
        aux_files += getInputFiles(outdir, "*.log", ignore_empty_files=False)
        aux_files += getInputFiles(".", "LikelihoodRatio.txt", ignore_empty_files=False)

        # Move the final groups file(s) to the groups dir
        groups_dir = makeDirOrdie("%s_groups_files" % outdir)
        bulk_move_to_dir(final_groups_files, groups_dir)

        # Move aux files to the aux dir
        aux_dir = makeAuxDir(outdir)
        bulk_move_to_dir(aux_files, aux_dir)

        # Cleanup the pool
        cleanup_pool(pool)


    def parseCROPoutToGroups(crop_out_file, output_groups_file):
        """Parses a CROP output file to a groups file.  Crop files are pretty close to the groups format, we just need to
            replace commas with spaces, and clip the counts from child names.
            e.g.
            Crop line:
                BALI4606_0_ID2033_1	BALI4606_0_ID2033_1,BALI4606_0_ID1668_1,BALI4606_0_ID2079_1

            groups line:
                BALI4606_0_ID2033   BALI4606_0_ID2033 BALI4606_0_ID1668 BALI4606_0_ID2079

        :param crop_out_file: The file path to the input crop output file.
        :param output_groups_file: The file path to write the output groups file to.
        :return: The filepath to the output groups file
        """

        with open(output_groups_file, 'w') as out:
            i = 0
            output = ""
            for line in open(crop_out_file, 'r'):
                i += 1
                if i % 100000 == 0:
                    out.write(output)
                    output = ""
                name, children = line.split("\t")
                seqs = [clip_count(seq) for seq in children.split(",")]
                output += "%s\t%s\n" % (clip_count(name), " ".join(seqs))
            out.write(output)
