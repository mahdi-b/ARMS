from classes.ChewbaccaProgram import *
from classes.Helpers import *
from classes.ProgramRunner import *

class Align_Program_Macse(ChewbaccaProgram):
    name = "macse"


    def execute_program(self):
        args = self.args
        self.align_macse(args.input_f, args.db, args.outdir, args.threads)


    def align_macse(self, input_f, db, outdir, threads):
        """Aligns sequences by iteratively adding them to a known good alignment.

        :param input_f: Filepath to an input file or folder to rename.
        :param db: Filepath to a reference file or folder of reference files for alignment.
        :param outdir: Filepath to the output directory.
        """
        # "macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
        #                                    \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
        #                                    -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
        #                                    -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",
        makeDirOrdie(outdir)
        inputs = getInputFiles(input_f)
        pool = init_pool(min(len(inputs), threads))
        printVerbose("Aligning reads using MACSE")
        inputs = getInputFiles(input_f)
        parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.MACSE_ALIGN,
                                                          [db, db, input] +
                                                          ["%s/%s" % (outdir, getFileName(input))] * 3,
                                                          {"exists": [input, db]}) for input in inputs], pool)
        printVerbose("Done with MACSE alignment.")
        cleanup_pool(pool)