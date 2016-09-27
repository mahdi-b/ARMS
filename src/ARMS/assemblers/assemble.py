from classes.Helpers import *
from classes.ProgramRunner import *


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
