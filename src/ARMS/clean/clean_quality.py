from classes.Helpers import *
from clean_quality_trimmomatic import clean_quality_trimmomatic

def clean_quality_main(input_f, outdir, program, threads, aux_params):
    """Removes areas of low quality in sequences from a FASTQ file.

    :param input_f: Filepath to input file or folder.
    :param outdir: Filepath to the output directory.
    :param program: Program to use to clean sequences.  Options are ["trimmomatic"].
    :param threads: Number of processes to use to clean the input fileset.
    :param aux_params: Program-specific named-parameters.
    """
    makeDirOrdie(outdir)
    inputs = getInputFiles(input_f)
    debugPrintInputInfo(inputs, "clean")
    pool = init_pool(min(len(inputs), threads))
    if program == "trimmomatic":
        clean_quality_trimmomatic(inputs, outdir, aux_params["window_size"], aux_params["quality"],
                                  aux_params["min_len"], pool)
    cleanup_pool(pool)



