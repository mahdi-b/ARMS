from classes.Helpers import *
from classes.ProgramRunner import *
from splitKperFasta import splitK


def partition_main(input_f, outdir, threads, program, chunksize, filetype, aux_params):
    """Partitions a fasta/fastq file into files with <chunksize> sequences.

    :param input_f: Filepath to a file or folder of files to partition.
    :param outdir: The directory to write split files to.
    :param threads: The number of processes to use to partition the input fileset.
    :param chunksize: The number of sequences per file.
    :param filetype: Either 'fasta' or 'fastq'.
    :param aux_params A dictionary of program-specific named-parameters.
    """
    makeDirOrdie(outdir)
    # Gather input files
    inputs = getInputFiles(input_f)
    debugPrintInputInfo(inputs, "partitioned")
    pool = init_pool(min(len(inputs), threads))
    if program == "chewbacca":
        partition_chewbacca(inputs, outdir, threads, chunksize, filetype, pool)
    cleanup_pool(pool)


def partition_chewbacca(inputs, outdir, threads, chunksize, filetype, pool):
    """Partition a fasta/fastq file into chunks of user-defined size.

    :param inputs: A list of filepaths to files to partition.
    :param outdir: The directory to write split files to.
    :param threads: The number of processes to use to partition the input fileset.
    :param chunksize: The number of sequences per file.
    :param filetype: Either 'fasta' or 'fastq'.
    :param pool: A fully-initalized multiprocessing.Pool object.
    """
    # def splitK(inputFasta, prefix, nbSeqsPerFile, filetype):
    printVerbose("Partitioning Files...")
    parallel(runPythonInstance,
             [(splitK, input_, "%s/%s" % (outdir, strip_ixes(input_)), chunksize, filetype)
              for input_ in inputs], pool)
    printVerbose("Done partitioning files.")