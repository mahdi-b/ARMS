from classes.ChewbaccaProgram import *
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from classes.Helpers import *
from cluster_helpers import handle_groups_file_update
from parse.parseCROPoutToGroups import parseCROPoutToGroups


class Cluster_Program_Crop(ChewbaccaProgram):
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

        :param input_f:
        :param outdir:
        :param groupsfile:
        :param threads:
        :param blocksize:
        :param clustpct:
        :param maxmcmc:
        :param maxsm:
        :param rare:
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