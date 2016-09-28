from classes.Helpers import *
from classes.ProgramRunner import *


def assemble_pear(inputs, outdir, name, pearthreads):
    """Uses PEAR to assemble paired F/R read files in parallel.

    :param inputs: A list of tuples of paried F/R sequence reads.
    :param outdir: File path to the output directory
    :param name: File prefix for the assembled reads files.
    :param pearthreads: The number of threads per process to use to assemble a pair of reads.
    """
    # "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s -m %d"
    makeDirOrdie(outdir)

    printVerbose("\tAssembling reads with pear")
    debugPrintInputInfo(inputs, "assemble")
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.ASSEMBLE_PEAR, [forwards, reverse,
                                                                                            "%s/%s_%s" % (
                                                                                                outdir, name,
                                                                                                getFileName(forwards)),
                                                                                            pearthreads],
                                                      {"exists": [forwards, reverse], "positive": [pearthreads]})
                                        for forwards, reverse in inputs], pool)

    printVerbose("Done assembling sequences...")
    # Grab all the auxillary files (everything not containing ".assembled."
    aux_files = getInputFiles(outdir, "*", "*.assembled.*", ignore_empty_files=False)
    # make aux dir for extraneous files and move them there
    bulk_move_to_dir(aux_files, makeAuxDir(outdir))