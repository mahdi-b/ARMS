from classes.Helpers import *
from assemble_pear import assemble_pear


def assemble_main(input_f, input_r, outdir, name, program, threads, aux_params):
    """Assembles reads from two (left and right) fastq files/directories.  For a set of k forward read files, and k
        reverse read files, return k assembled files.  Matching forward and reverse files should be identically named,
        except for a <forward>/<reverse> suffix that indicates the read orientation.  Two suffix pairs are supported:
        '_forwards' and '_reverse',
        and
        '_R1' and 'R2'
        Choose ONE suffix style and stick to it.
        e.g. Sample_100_forwards.fq and Sample_100_reverse.fq will be assembled into Sample_100_assembled.fq.
          Alternatively, Sample_100_R1.fq and Sample_100_R2.fq will be assembled into Sample_100_assembled.fq.
          You can provide as many pairs of files as you wish as long as they follow exactly on of the above naming
          conventions.  If a 'name' parameter is provided, it will be used as a suffix for all assembled sequence files.

    :param input_f: File path to forward Fastq Reads file or folder.
    :param input_r: File path to reverse Fastq Reads file or folder.
    :param outdir: File path to the output directory.
    :param name: File prefix for the assembled reads.
    :param program: The program to use to complete the assembly.  Choices are ["pear"].
    :param threads: The number of processes to use to complete assembly.
    :param aux_params: Program-specific named-parameters.
    """
    # Collect input files, and validate that they exist and match
    inputs = validate_paired_fastq_reads(input_f, input_r)
    pool = init_pool(min(len(inputs), threads))
    if program == "pear":
        pearthreads = int(aux_params['pearthreads'])
        assemble_pear(inputs, outdir, name, threads, pearthreads)
    cleanup_pool(pool)



