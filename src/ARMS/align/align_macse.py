from classes.Helpers import *
from classes.ProgramRunner import *

def align_macse(input_f, ref, outdir, pool):
    """Aligns sequences by iteratively adding them to a known good alignment.

    :param input_f: Filepath to an input file or folder to rename.
    :param ref: Filepath to a reference file or folder of reference files for alignment.
    :param outdir: Filepath to the output directory.
    :param pool: A fully initalized multiprocessing.Pool object.
    """
    # "macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
    #                                    \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
    #                                    -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
    #                                    -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",
    printVerbose("Aligning reads using MACSE")
    inputs = getInputFiles(input_f)
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.MACSE_ALIGN,
                                                      [ref, ref, input] +
                                                      ["%s/%s" % (outdir, getFileName(input))] * 3,
                                                      {"exists": [input, ref]}) for input in inputs], pool)
    printVerbose("Done with MACSE alignment.")