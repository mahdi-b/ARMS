from align.Align_Command import align_main
from align.align_clean import align_clean_main
from assemble.Assemble_Command import Assemble_Command
from classes.Helpers import *
from clean.Clean_Adapters_Command import Clean_Adapters_Command
from clean.Clean_Quality_Command import Clean_Quality_Command
from clean.clean_gap_chars import ungap_main
from cluster.cluster import cluster_main
from demux.Demux_Command import Demux_Command
from dereplicate.dereplicate import dereplicate_main
from preclean.Preclean_Command import *
from merge.merge import merge_main
from otu.annotate import annotate_main
from otu.build import build_otu_main
from otu.query import query_main, QUERY_TYPE
from rename.Rename_Command import Rename_Command
from util.fastqToFasta import translateFastqToFasta
from util.partition import partition_main


# TODO: reconcile args parameter documentation
def preclean(args):
    """Attempts to fix errors in reads before assembly.
    :param args: An argparse object with the following parameters:
                    name            Textual ID for the data set.
                    input_f         Forward Fastq Reads file or directory.
                    input_r         Reverse Fastq Reads file or directory.
                    threads         The number of threads to use durring assembly.
                    outdir          Directory where outputs will be saved.
                    program         The program to use for precleaning.  Choices are: ['bayeshammer']
                                    Default is bayeshammer.
    """
    Preclean_Command(args).execute_command()


def assemble(args):
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
    :param args: An argparse object with the following parameters:
                    name            Textual ID for the data set.
                    input_f         Forward Fastq Reads file or directory.
                    input_r         Reverse Fastq Reads file or directory.
                    threads         The number of threads to use durring assembly.
                    outdir          Directory where outputs will be saved.
    """

    Assemble_Command(args).execute_command()


def demultiplex(args):
    """Splits a fasta/fastq file on a set of barcodes.  For a set of k input files, each input file is assigned a
    file_id. Then, each file is split into a set of subfiles, where each subfile named <sample>_splitOut_<k> contains
    the sequences belonging to <sample> from file <k>.
    Note, the assignment of a file to its file_id is arbitrary and should not be used to identify file groupings.
    E.x:
    File Alpha.fastq contains sequences belonging to samples A,B, and C,
    File Beta.fastq contains sequences belonging to samples A and C,
    we expect to see the output files:
    A_splitOut_0.fastq #sequences from sample A that were in Alpha
    A_splitOut_1.fastq #sequences from sample A that were in Beta
    B_splitOut_0.fastq #sequences from sample B that were in Alpha
    C_splitOut_0.fastq #sequences from sample C that were in Alpha
    C_splitOut_1.fastq #sequences from sample C that were in Beta

    :param args: An argparse object with the following parameters:
                    input       Input file or folder containing input files.
                    barcodes    Tab delimited file mapping barcodes to their samples.  Must be a single file.
                    outdir      Directory where outputs will be saved.
    """
    Demux_Command(args).execute_command()


def rename_sequences(args):
    """Renames sequences in a fasta/fastq file as <filename>_ID0, <filename>_ID1, <filename>_ID2, etc., where <filename>
        is the name of the fasta/fastq file without any extensions or chewbacca suffixes.
    :param args: An argparse object with the following parameters:
                    input       Input file or folder containing only fasta or fastq files.
                    outdir      Directory where outputs will be saved.
                    filetype    File type of the input files.  Either 'fasta' or 'fastq'.
                    clip        True if input file names contain trailing demux_seqs identifiers.
                                    e.g. True if file name were 'Sample_395_0.split.out', 'Sample_395_1.split.out', etc.
                                    Important because sample names are derived from input file names.
    """
    Rename_Command(args).execute_command()


def trim_adapters(args):
    """Trim adapters (and preceeding barcodes) from sequences in input file(s).  Sequences should be in
        the following format: <BARCODE><ADAPTER><SEQUENCE><RC_ADAPTER>, where ADAPTER is defined in the adapters file,
        and RC_ADAPTER is defined in the rcadapters file.  By default, Flexbar does not allow any 'N' characters in
        SEQUENCE, and will toss any sequences that do contain 'N'.  To avoid this, use the -u or --allowedns flags to
        specify the maximum number of 'N's to allow
    :param args: An argparse object with the following parameters:
                    input       Input file or folder containing only fasta or fastq files.
                    outdir      Directory where outputs will be saved.
                    adapters    Filepath to a single file listing forwards adapters, or directory of such file(s).
                    adaptersrc  Filepath to a single file listing reverse complemented adapters, or directory of such
                                    file(s).
                    allowedns   Non-negative integer value indicating the maximum number of 'N's to tolerate in a
                                    sequence.
    """
    Clean_Adapters_Command(args).execute_command()


def clean_quality(args):
    """Uses a sliding window to identify and trim away areas of low quality.

    :param args: An argparse object with the following parameters:
                    inputFile	Input Fastq file
                    outputFile	Output Fastq file
                    windowSize	Width of the sliding window
                    quality 	Minimum passing quality for the sliding window
                    minLen	    Minimum passing length for a cleaned sequence
    """
    Clean_Quality_Command(args).execute_command()


def dereplicate(args):
    """Dereplicates a set of fasta files.  Identical (or identical spanning) sequences are considered
        replicants.  (100% match).  NOTE: only dereplicates within each fasta file (not across all files).  Merge
        files before hand if you want to dereplciate across multiple files.

    :param args: An argparse object with the following parameters:
                    input_f:        Filepath to the file or folder of files to dereplicate.
                    outdir:         Filepath to the output directory.
                    groupsfile:     A groups file to use as a reference for dereplication counting.  If no groups file
                                            is provided, input sequences are conidered singletons (regardless of their
                                            dereplication count).
                    threads:        The number of processes to use to dereplicate the fileset.
                    stripcounts:    If True, strips the trailing dereplication counts from a file before dereplication.
                    aux_args:       A dictionary of program-specific named-parameters.
    """
    dereplicate_main(args.input_f, args.outdir, args.groupsfile, args.program, args.threads, args.stripcounts,
                     vars(args))


def partition(args):
    """ Partition a fasta/fastq file into multiple files, each with <chunksize> sequences.

    :param args: An argparse object with the following parameters:
                    input	    Input fasta file to split
                    outdir      Directory where outputs will be saved
                    chunksize	Chunksize.
                    filetype	Filetype of the files to be partitioned
    """
    partition_main(args.input_f, args.outdir, args.threads, args.program, args.chunksize, args.filetype, vars(args))


def merge(args):
    """Concatenates files in a directory.
       :param args: An argparse object with the following parameters:
                       input        Cleaned inputs File.
                       outdir       Directory where outputs will be saved.
                       name         Name prefix for the merged file.
                       fileext      Output file extension.  e.g 'fasta', 'fastq', 'txt'
       """
    merge_main(args.input_f, args.outdir, args.program, args.name, args.fileext, vars(args))


def ungap_fasta(args):
    """Removes a gap character from sequences in a fasta file.  Useful for removing characters from an alignment file.
       :param args: An argparse object with the following parameters:
                       input        Fasta files to ungap.
                       outdir       Directory where outputs will be saved.
                       gapchar      A gap character to remove from the sequences.
                       fileext      Either 'fastq' or 'fasta'.
    """
    ungap_main(args.input_f, args.outdir, args.gapchars, args.fileext, args.program, args.threads, vars(args))


def cluster(args):
    """Clusters a fasta file.

    :param args: An argparse object.
    """
    cluster_main(args.input_f, args.outdir, args.program, args.groupsfile, args.threads, vars(args))


def build_matrix(args):
    """Builds the unannotated OTU table.
    :param args: An argparse object with the following parameters:
                    groups       File/folder with .groups files.
                    outdir      Directory to put the output files.
                    samples      File/folder with .samples files.
                    barcodes    A single barcodes file listing all possible sample groups.
    """

    build_otu_main(args.outdir, args.groups_file, args.samples_file, args.program, args.barcodes_file, vars(args))


def query_fasta(args):
    """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

    :param args: An argparse object with the following parameters:
                    accnosFile  List of sequence names to remove
                    outdir      Directory to put the output files
    """
    query_main(args.input_f, args.outdir, args.threads, args.program, QUERY_TYPE.FASTA, vars(args))


def query_ncbi(args):
    """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

    :param args: An argparse object with the following parameters:
                    input       Input file/folder with fasta sequences
                    outdir      Directory to put the output files
    """
    query_main(args.input_f, args.outdir, args.threads, args.program, QUERY_TYPE.DATABASE, vars(args))


def annotate_matrix(args):
    """Annotates an OTU table.
    :param args: An argparse object with the following parameters:
                    input       File/folder with matrix files.
                    outdir      Directory to put the output files.
                    map      File/folder with .samples files.
    """
    annotate_main(args.input_f, args.annotation, args.outdir, args.threads, args.program, vars(args))


def make_fasta(args):
    """Converts a fastq file to fasta format.
    :param args: An argparse object with the following parameters:
                    input       File/folder with fastq files.
                    outdir      Directory to put the output files.xs
    """
    makeDirOrdie(args.outdir)
    inputs = getInputFiles(args.input_f)
    debugPrintInputInfo(inputs, "convert to fasta.")
    printVerbose("Converting to fasta...")

    pool = init_pool(min(len(inputs), args.threads))
    parallel(runPythonInstance, [(translateFastqToFasta, input_, "%s/%s.fasta" % (args.outdir, getFileName(input_)))
                                 for input_ in inputs], pool)
    printVerbose("Done converting.")
    cleanup_pool(pool)


def mase_align(args, pool=Pool(processes=1)):
    """Aligns sequences by iteratively adding them to a known good alignment.
     :param args: An argparse object with the following parameters:
                    db                  Database against which to align and filter reads
                    samplesDir          Directory containig the samples to be cleaned
                    outdir              Directory where outputs will be saved
     :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    align_main(args.input_f, args.ref, args.outdir, args.program, args.threads, vars(args))


def clean_macse_align(args, pool=Pool(processes=1)):
    """Removes non-nucleotide characters in MACSE aligned sequences for all fasta files in the samples directory
        (the samplesDir argument).
    :param args: An argparse object with the following parameters:
                    samplesDir          Directory containig the samples to be cleaned
                    outdir              Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    align_clean_main(args.input_f, args.samplesdir, args.ref, args.outdir, args.threads, args.program, vars(args))
