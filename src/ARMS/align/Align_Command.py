from classes.Helpers import *
from align_macse import align_macse


# TODO: While vsearch is an alignment tool, we never use it as such, and don't have code snippets to support it as a
# TODO:  multi-sequence aligner.  Therefore, its included as part of the "query" step.  Look for it in otu/query.py
def align_main(input_f, ref, outdir, program, threads, aux_params):
    """Aligns input sequence file(s) against known good sequences.

    :param input_f: Filepath to an input file or folder to rename.
    :param ref: Filepath to a reference file or folder of reference files for alignment.
    :param outdir: Filepath to the output directory.
    :param program: The program to use to align sequences.  Choices: ['macse'].  Default: 'macse'.
    :param threads: Number of processes to use to align the input fileset.
    :param aux_params: A dictionary of porgram-specific named-parameters.
    """
    makeDirOrdie(outdir)
    inputs = getInputFiles(input_f)
    pool = init_pool(min(len(inputs), threads))
    if program == "macse":
        align_macse(input_f, ref, outdir, pool)
    cleanup_pool(pool)


