from classes.Helpers import *
from classes.ProgramRunner import *
from align_clean_macse import align_clean_macse


def align_clean_main(input_f, samplesdir, ref, outdir, threads, program, aux_params):
    """
    :param input_f: File path to file or folder of files to clean.
    :param samplesdir: Filepath to the original, unaligned input files (the inputs to the macse aligner).
    :param ref: Filepath to the reference file used to align the input files.
    :param outdir: Filepath to the directory to write outputs to.
    :param threads: The number of processes to use to clean the alignments.
    :param program: The program to use to clean the alignments.  Choices are ["macse"].  Default is "macse".
    :param aux_params: A dictionary of program-specific named-parameters.
    """
    makeDirOrdie(outdir)
    inputs = getInputFiles(input_f)
    pool = init_pool(min(len(inputs), threads))
    if program == "macse":
        align_clean_macse(input_f, ref, samplesdir, outdir, pool)
    cleanup_pool(pool)
