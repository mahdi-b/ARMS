from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import getInputFiles, init_pool, printVerbose, run_parallel, getFileName, cleanup_pool
from classes.ProgramRunner import ProgramRunner,ProgramRunnerCommands

class Clean_Deep_Program_Macse(ChewbaccaProgram):
    """Uses Macse's enrichAlignment program to align a set of sequences."""
    name = "macse"

    def execute_program(self):
        args = self.args
        self.align_macse(args.input_f, args.db, args.outdir, args.processes, args.extraargstring)

    def align_macse(self, input_f, db, outdir, processes, extraargstring):
        """Aligns sequences by iteratively adding them to a known good alignment.

        :param input_f: Filepath to an input file or folder to rename.
        :param db: Filepath to a reference file or folder of reference files for alignment.
        :param outdir: Filepath to the output directory.
        :param processes: The maximum number of processes to use.
        :param extraargstring: Advanced program parameter string.
        """
        # "macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
        #                                    \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
        #                                    -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
        #                                    -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",

        inputs = getInputFiles(input_f)
        pool = init_pool(min(len(inputs), processes))
        printVerbose("Aligning reads using MACSE")
        inputs = getInputFiles(input_f)
        run_parallel([ProgramRunner(ProgramRunnerCommands.MACSE_ALIGN,
                                    [db, db, input_] + ["%s/%s" % (outdir, getFileName(input_))] * 3,
                                    {"exists": [input_, db]}, extraargstring)
                      for input_ in inputs], pool)
        printVerbose("Done with MACSE alignment.")
        cleanup_pool(pool)
