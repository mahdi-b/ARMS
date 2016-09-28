from classes.Helpers import *
from dereplicate_vsearch import dereplicate_vsearch


def dereplicate_main(input_f, outdir, groupsfile, program, threads, stripcounts, aux_args):
    """Dereplicates a set of fasta files.  Identical (or identical spanning) sequences are considered \
        replicants.  (100% match).  NOTE: only dereplicates within each fasta file (not across all files).  Merge \
        files before hand if you want to dereplciate across multiple files.

    :param input_f: Filepath to the file or folder of files to dereplicate.
    :param outdir: Filepath to the output directory.
    :param groupsfile: A groups file to use as a reference for dereplication counting.  If no groups file is provided,
                        input sequences are conidered singletons (regardless of their dereplication count).
    :param threads: The number of processes to use to dereplicate the fileset.
    :param stripcounts: If True, strips the trailing dereplication counts from a file before dereplication.
    :param aux_args: A dictionary of program-specific named-parameters.
    """
    makeDirOrdie(outdir)
    inputs = getInputFiles(input_f)
    pool = init_pool(min(len(inputs), threads))

    if program == "vsearch":
        dereplicate_vsearch(inputs, outdir, groupsfile, threads, stripcounts, pool)
    cleanup_pool(pool)


