from classes.Helpers import *
from preclean_bayeshammer import preclean_bayeshammer

def preclean_main(input_f, input_r, outdir, program, threads, aux_params):
    """Assembles reads from two (left and right) fastq files/directories.  Switchboard function to call the appropriate
        subprogram.

    :param input_f: File path to file or folder of left reads to clean.
    :param input_r: File path to file or folder of right reads to clean.
    :param outdir: Filepath to output directory.
    :param program: The program to use.  Current options are "bayeshammer".
    :param threads: Number of processes to use to process the input files.
    :param aux_params: dictionary of program-specific commands.
    """

    makeDirOrdie(outdir)
    # Collect input files, and validate that they match
    inputs = validate_paired_fastq_reads(input_f, input_r)
    pool = init_pool(min(len(inputs), threads))
    if program == "bayeshammer":
        return preclean_bayeshammer(inputs, outdir, threads, pool, aux_params["bayesthreads"])
    cleanup_pool(pool)
