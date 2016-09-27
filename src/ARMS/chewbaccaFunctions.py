from itertools import product

from assemblers.assemble import assemble_main
from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from cluster.cluster import cluster_main
from converters.clean_macse_alignment import clean_macse_alignment
from converters.fastqToFasta import translateFastqToFasta
from converters.ungap import ungap
from demultiplexers.demultiplex import demultiplex_main
from otu_tables.annotateOTUtable import annotateOTUtable
from otu_tables.buildOTUtable import buildOTUtable
from parsers.parseUCtoGroups import parseUCtoGroups
from parsers.parseVSearchoutForTaxa import *
from renamers.rename import rename_main
from renamers.renameWithReplicantCounts import renameWithReplicantCounts
from renamers.renameWithoutCount import removeCountsFromFastFile
from renamers.updateGroups import update_groups
from src.ARMS.fixers.error_fix import preclean_main
from trimmers.trimmers import trim_adapters_main
from utils.joinFiles import joinFiles
from utils.splitKperFasta import splitK


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
    preclean_main(args.input_f, args.input_r, args.outdir, args.program, args.threads, vars(args))


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

    assemble_main(args.input_f, args.input_r, args.outdir, args.name, args.program, args.threads, vars(args))


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
    # "barcode.splitter": "cat \"%s\" | " + programPaths["FASTX"] + 'fastx_barcode_splitter.pl  --bcfile "%s" \
    #                                    -prefix "%s" --suffix %s --bol --mismatches 1',
    demultiplex_main(args.input_f, args.barcodes, args.outdir, args.programs, args.threads, vars(args))


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

    rename_main(args.input_f, args.outdir, args.filetype, args.clip, args.program, args.threads, vars(args))


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
    trim_adapters_main(args.input_f, args.adapters, args.adaptersrc, args.outdir, args.program, args.threads,
                       vars(args))


def trimmomatic(args):
    """Uses a sliding window to identify and trim away areas of low quality.

    :param args: An argparse object with the following parameters:
                    inputFile	Input Fastq file
                    outputFile	Output Fastq file
                    windowSize	Width of the sliding window
                    quality 	Minimum passing quality for the sliding window
                    minLen	    Minimum passing length for a cleaned sequence
    """
    # "trimomatic":       "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
    # -%phred %input %output SLIDINGWINDOW:%windowsize:%minAvgQuality MINLEN:%minLen"
    makeDirOrdie(args.outdir)

    printVerbose("Cleaning sequences with Trimmomatic...")
    inputs = getInputFiles(args.input_f)
    pool = init_pool(min(len(inputs), args.threads))

    debugPrintInputInfo(inputs, "clean")
    parallel(runProgramRunnerInstance,
             [ProgramRunner(ProgramRunnerCommands.CLEAN_TRIMMOMATIC,
                            [input_, "%s/%s_cleaned.fastq" % (args.outdir, strip_ixes(input_)),
                             args.windowSize, args.quality,
                             args.minLen],
                            {"exists": [args.outdir, input_],
                             "positive": [args.windowSize, args.quality, args.minLen]})
              for input_ in inputs], pool)
    printVerbose("Done cleaning sequences.")

    cleanup_pool(pool)


def dereplicate(args):
    makeDirOrdie(args.outdir)
    inputs = getInputFiles(args.input_f)
    pool = init_pool(min(len(inputs), args.threads))
    # REMOVES COUNTS FROM SEQUENCE NAMES IN ORDER TO CLUSTER PROPERLY
    # strip counts if we need to.
    if args.stripcounts:
        printVerbose("Removing counts from sequence names...")
        debugPrintInputInfo(inputs, "renamed")
        parallel(runPythonInstance, [(removeCountsFromFastFile, input_,
                                      "%s/%s_uncount.fasta" % (args.outdir, strip_ixes(input_)), 'fasta')
                                     for input_ in inputs], pool)
        printVerbose("Done removing counts.")

        # Grab the cleaned files as input for the next step
        inputs = getInputFiles(args.outdir, "*_uncount.fasta")

    # DEREPLICATE ONE MORE TIME
    printVerbose("Dereplicating before clustering...")
    debugPrintInputInfo(inputs, "dereplicated")
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.DEREP_VSEARCH,
                                                      [args.threads, input_,
                                                       "%s/%s_derep.fasta" % (args.outdir, strip_ixes(input_)),
                                                       "%s/%s_uc.out" % (args.outdir, strip_ixes(input_))],
                                                      {"exists": [input_], "positive": [args.threads]})
                                        for input_ in inputs], pool)
    printVerbose("Done dereplicating")

    # LOG DEREPLICATED SEQUENCES INTO A .GROUPS FILE
    # generates a .groups file named _uc_parsed.out
    # python parseUCtoGroups.py uc.out uc_parsed.out
    input_ucs = getInputFiles(args.outdir, "*_uc.out")
    printVerbose("Generating a groups file from dereplication.")
    debugPrintInputInfo(inputs, "parsed (into agroups file)")
    parallel(runPythonInstance,
             [(parseUCtoGroups, input_, "%s/%s_derep.groups" % (args.outdir, strip_ixes(input_)))
              for input_ in input_ucs], pool)

    most_recent_groups_files = getInputFiles(args.outdir, "*_derep.groups")

    # UPDATE THE MOST CURRENT GROUPS FILES WITH DEREPLICATION COUNTS
    if args.groupsfile is not None:
        # Grab the oldgroups file and the dereplicated groups file
        old_groups_files = getInputFiles(args.groupsfile)
        derep_groups_files = getInputFiles(args.outdir, "*_derep.groups")

        printVerbose("Updating .groups files with dereplicated data")
        logging.debug("%d Reference (old)groups files to be read:" % len(old_groups_files))
        logging.debug(str(old_groups_files))
        logging.debug("%d Dereplicated (new)groups files to be read:" % len(derep_groups_files))
        logging.debug(str(derep_groups_files))
        # update_groups (old_groups_files, new_groups_files, updated)
        update_groups(old_groups_files, derep_groups_files, args.outdir, "dereplicated")
        most_recent_groups_files = getInputFiles(args.outdir, "dereplicated*")
        printVerbose("Done updating .groups files.")

    if len(inputs) != len(most_recent_groups_files):
        print ("Error: Number of input fastas (%d) is not equal to the number ofgroups files (%d)." %
               (len(inputs), len(most_recent_groups_files)))
        exit()
    fasta_groups_pairs = zip(inputs, most_recent_groups_files)
    # ADD COUNT TO SEQUENCE NAMES AND SORT BY COUNT
    # python renameWithReplicantCounts.py
    #               8_macse_out/MACSEOUT_MERGED.fasta uc_parsed.out dereplicated_renamed.fasta
    printVerbose("Adding dereplication data to unique fasta")
    parallel(runPythonInstance,
             [(renameWithReplicantCounts, fasta, groups,
               "%s/%s_counts.fasta" % (args.outdir, strip_ixes(fasta)), 'fasta')
              for fasta, groups in fasta_groups_pairs], pool)
    printVerbose("Done adding data")

    aux_dir = makeAuxDir(args.outdir)
    groups_dir = makeDirOrdie("%s_groups_files" % args.outdir)
    bulk_move_to_dir(most_recent_groups_files, groups_dir)
    aux_files = getInputFiles(args.outdir, '*', "*_counts.fasta")
    bulk_move_to_dir(aux_files, aux_dir)


def partition(args):
    """ Partition the cleaned file into chunks of user-defined size.

    :param args: An argparse object with the following parameters:
                    input	    Input fasta file to split
                    outdir      Directory where outputs will be saved
                    chunksize	Chunksize.
                    filetype	Filetype of the files to be partitioned
    """
    # def splitK(inputFasta, prefix, nbSeqsPerFile, filetype):
    makeDirOrdie(args.outdir)

    # Gather input files
    inputs = getInputFiles(args.input_f)
    debugPrintInputInfo(inputs, "partitioned")
    printVerbose("Partitioning Files...")

    pool = init_pool(min(len(inputs), args.threads))

    parallel(runPythonInstance,
             [(splitK, input_, "%s/%s" % (args.outdir, strip_ixes(input_)), args.chunksize, args.filetype)
                for input_ in inputs], pool)
    printVerbose("Done partitioning files.")
    cleanup_pool(pool)


def merge(args):
    """Blindly concatenates files in a directory.
       :param args: An argparse object with the following parameters:
                       input        Cleaned inputs File.
                       outdir       Directory where outputs will be saved.
                       name         Name prefix for the merged file.
                       fileext      Output file extension.  e.g 'fasta', 'fastq', 'txt'
       """
    makeDirOrdie(args.outdir)
    inputs = getInputFiles(args.input_f)
    debugPrintInputInfo(inputs, "merged together.")
    printVerbose("Merging files...")
    joinFiles(inputs, "%s/%s_MERGED.%s" % (args.outdir, args.name, args.fileext))
    pool = init_pool(min(len(inputs), args.threads))
    printVerbose("Done merging.")
    cleanup_pool(pool)


def ungap_fasta(args):
    """Removes a gap character from sequences in a fasta file.  Useful for removing characters from an alignment file.
       :param args: An argparse object with the following parameters:
                       input        Fasta files to ungap.
                       outdir       Directory where outputs will be saved.
                       gapchar      A gap character to remove from the sequences.
                       fileext      Either 'fastq' or 'fasta'.
       """
    makeDirOrdie(args.outdir)
    input_files = getInputFiles(args.input_f,"*.fasta")
    debugPrintInputInfo(input_files, "ungap.")
    pool = init_pool(min(len(input_files), args.threads))
    printVerbose("Removing all '%s' from sequences..." % args.gapchar)
    # ungap(file_to_clean, output_file_name, gap_char, file_type):
    parallel(runPythonInstance, [(ungap, input_, "%s/%s_cleaned.%s" % (args.outdir, strip_ixes(input_), 'fasta'),
                                  args.gapchar, args.fileext) for input_ in input_files], pool)
    printVerbose("Done removing.")
    cleanup_pool(pool)


def cluster(args):
    """Clusters a fasta file.

    :param args: An argparse object.
    """
    cluster_main(args.input_f, args.outdir, args.program, args.groupsfile, args.threads, vars(args))


# TODO doc and debug
def build_matrix(args):
    """Builds the unannotated OTU table.
    :param args: An argparse object with the following parameters:
                    groups       File/folder with .groups files.
                    outdir      Directory to put the output files.
                    samples      File/folder with .samples files.
                    barcodes    A single barcodes file listing all possible sample groups.
    """
    makeDirOrdie(args.outdir)
    groups = getInputFiles(args.groups)
    debugPrintInputInfo(groups, "read.")
    samples = getInputFiles(args.samples)
    debugPrintInputInfo(samples, "read.")
    barcodes = getInputFiles(args.barcodes)
    debugPrintInputInfo(barcodes, "read.")
    printVerbose("Building matrix...")
    buildOTUtable(groups, samples, barcodes[0], "%s/%s.txt" % (args.outdir, "matrix"))
    printVerbose("Done building.")


def query_fasta(args):
    """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

    :param args: An argparse object with the following parameters:
                    accnosFile  List of sequence names to remove
                    outdir      Directory to put the output files
    """
    aln_user_string = ""
    makeDirOrdie(args.outdir)

    # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
    #       --userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt

    # expecting a fasta to annotate
    query_fastas = getInputFiles(args.input_f)
    debugPrintInputInfo(query_fastas, "queried for identification.")
    ref_fastas = getInputFiles(args.referencefasta)
    debugPrintInputInfo(ref_fastas, "referenced for sequence identification.")
    tax_info_files = getInputFiles(args.taxinfo)
    debugPrintInputInfo(tax_info_files, "referenced for taxanomic names.")

    # make sure the number of reference fasta files is the same as the number of tax_info files
    if len(tax_info_files) != len(ref_fastas):
        print "Error: The number of reference fastas and taxonomic mapping files is not the same.  There must be one \
                taxonomic mapping file for each reference fasta."
        sys.exit()

    ref_data_pairs = zip(ref_fastas, tax_info_files)
    inputs = [x for x in product(query_fastas, ref_data_pairs)]

    pool = init_pool(min(len(query_fastas), args.threads))
    printVerbose("Querying...")
    # input, db, userout, alnout, aln_userfield
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.ALIGN_VSEARCH,
                                                      [args.threads, query, ref_fasta,
                                                       "%s/%s.out" % (args.outdir, strip_ixes(query)),
                                                       "%s/%s.alnout" % (args.outdir, strip_ixes(query)),
                                                       aln_user_string],
                                                      {"exists": [query, ref_fasta], "positive": [args.threads]})
                                        for query, (ref_fasta, tax_info) in inputs], pool)
    printVerbose("Done querying.")
    printVerbose("Parsing output...")

    # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
    # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
    #
    # parseVSearchOutputAgainstFasta(vsearch_outfile, taxInfo, output_file, min_simmilarity, min_coverage):
    parallel(runPythonInstance, [(parseVSearchOutputAgainstFasta, "%s/%s.out" % (args.outdir, strip_ixes(query)),
                                  tax_info, "%s/%s_result.out" % (args.outdir, strip_ixes(query)),
                                  args.simmilarity, args.coverage)
                                 for query, (ref_fasta, tax_info) in inputs], pool)
    printVerbose("\nDone parsing...")
    # Gather and move auxillary files
    aux_files = getInputFiles(args.outdir, "*", "*_result.out")
    bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    cleanup_pool(pool)


def query_ncbi(args):
    """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

    :param args: An argparse object with the following parameters:
                    input       Input file/folder with fasta sequences
                    outdir      Directory to put the output files
    """
    coi_fasta = os.path.expanduser("~/ARMS/refs/COI.fasta")
    ncbi_db_string = os.path.expanduser("~/ARMS/refs/ncbi.db")
    aln_user_string = "--userfields query+target+id+alnlen+qcov"
    makeDirOrdie(args.outdir)

    # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
    # --userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt

    # expecting a fasta to annotate
    query_fastas = getInputFiles(args.input_f)
    debugPrintInputInfo(query_fastas, "queried agains NCBI.")
    pool = init_pool(min(len(query_fastas), args.threads))
    printVerbose("Querying NCBI...")
    # input, db, userout, alnout, aln_userfield
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.ALIGN_VSEARCH,
                                                      [args.threads, query, coi_fasta,
                                                       "%s/%s.out" % (args.outdir, strip_ixes(query)),
                                                       "%s/%s.alnout" % (args.outdir, strip_ixes(query)),
                                                       aln_user_string],
                                                      {"exists": [query, coi_fasta], "positive": [args.threads]})
                                        for query in query_fastas], pool)

    printVerbose("Done with query.")
    printVerbose("Parsing output...")
    # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
    # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
    #
    # parseVSearchOutputAgainstNCBI(vsearch_out, ncbi_db, min_coverage, min_similarity)> parsed_nt.out
    parallel(runPythonInstance,
             [(parseVSearchOutputAgainstNCBI, "%s/%s.out" % (args.outdir, strip_ixes(query)), ncbi_db_string,
               "%s/%s_result.out" % (args.outdir, strip_ixes(query)),
               97, 85) for query in query_fastas], pool)
    printVerbose("Done processing.")

    # Gather and move auxillary files
    aux_dir = makeAuxDir(args.outdir)
    aux_files = getInputFiles(args.outdir, "*", "*_result.out")
    bulk_move_to_dir(aux_files, aux_dir)

    cleanup_pool(pool)


# TODO doc and debug
def annotate_matrix(args):
    """Annotates an OTU table.
    :param args: An argparse object with the following parameters:
                    input       File/folder with matrix files.
                    outdir      Directory to put the output files.
                    map      File/folder with .samples files.
    """
    makeDirOrdie(args.outdir)
    matricies = getInputFiles(args.input_f)
    debugPrintInputInfo(matricies, "annotated.")
    annotations = getInputFiles(args.annotation)
    debugPrintInputInfo(annotations, "parse.")
    inputs = product(matricies, annotations)

    pool = init_pool(min(len(matricies) * len(annotations), args.threads))
    printVerbose("Annotating matrix...")
    parallel(runPythonInstance, [(annotateOTUtable, matrix, annotation, "%s/%s.txt" % (args.outdir, "matrix"))
                                 for matrix, annotation in inputs], pool)
    printVerbose("Done Annotating.")

    cleanup_pool(pool)


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

def macseAlignSeqs(args, pool=Pool(processes=1)):
    """Aligns sequences by iteratively adding them to a known good alignment.
     :param args: An argparse object with the following parameters:
                    db                  Database against which to align and filter reads
                    samplesDir          Directory containig the samples to be cleaned
                    outdir              Directory where outputs will be saved
     :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # "macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
    #                                    \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
    #                                    -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
    #                                    -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",
    makeDirOrdie(args.outdir)

    printVerbose("\t %s Aligning reads using MACSE")
    inputs = getInputFiles(args.input_f)
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.MACSE_ALIGN,
                                              [args.db, args.db, input] +
                                              ["%s/%s" % (args.outdir, getFileName(input))] * 3,
                                              {"exists": [input, args.db]}) for input in inputs], pool)
    cleanup_pool(pool)


def macseCleanAlignments(args, pool=Pool(processes=1)):
    """Removes non-nucleotide characters in MACSE aligned sequences for all fasta files in the samples directory
        (the samplesDir argument).
    :param args: An argparse object with the following parameters:
                    samplesDir          Directory containig the samples to be cleaned
                    outdir              Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # "macse_format":     "java -jar " + programPaths["MACSE"] + "  -prog exportAlignment -align \"%s\" \
    #
    #                                  -charForRemainingFS - -gc_def 5 -out_AA \"%s\" -out_NT \"%s\" -statFile \"%s\""
    makeDirOrdie(args.outdir)
    printVerbose("\t %s Processing MACSE alignments")
    samplesList = getInputFiles(args.samplesdir)
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.MACSE_FORMAT,
                                              ["%s/%s_NT" % (args.input_f, strip_ixes(sample)),
                                               "%s/%s_AA_macse.fasta" % (args.outdir, strip_ixes(sample)),
                                               "%s/%s_NT_macse.fasta" % (args.outdir, strip_ixes(sample)),
                                               "%s/%s_macse.csv" % (args.outdir, strip_ixes(sample))],
                                              {"exists": []}) for sample in samplesList], pool)
    printVerbose("\tCleaning MACSE alignments")

    printVerbose("Processing %s samples..." % len(samplesList))
    nt_macse_outs = ["%s/%s_NT_macse.fasta" % (args.outdir, strip_ixes(sample)) for sample in samplesList]

    # Clean the alignments
    parallel(runPythonInstance, [(clean_macse_alignment, input_, args.db,
                                  "%s/%s" % (args.outdir,"%s_cleaned.fasta" % strip_ixes(input_)))
                                 for input_ in nt_macse_outs], pool)

    # Cat the cleaned alignments
    cleaned_alignments = getInputFiles(args.outdir, "*_cleaned.fasta")
    joinFiles(cleaned_alignments,"%s/MACSE_OUT_MERGED.fasta" % args.outdir)

    aux_dir = makeAuxDir(args.outdir)
    aux_files = getInputFiles(args.outdir, "*", "MACSE_OUT_MERGED.fasta")
    bulk_move_to_dir(aux_files, aux_dir)
    cleanup_pool(pool)