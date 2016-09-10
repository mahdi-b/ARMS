from itertools import product
from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from converters.capitalizeSeqs import capitalizeSeqs
from converters.fastqToFasta import translateFastqToFasta
from converters.seedToNames import seedToNames
from converters.ungap import ungap
from otu_tables.annotateOTUtable import annotateOTUtable
from otu_tables.buildOTUtable import buildOTUtable
from parsers.parseUCtoNames import parseUCtoNames
from parsers.parseVSearchoutForTaxa import parseVSearchOutputAgainstFasta, parseVSearchOutputAgainstNCBI
from renamers.renameSerially import serialRename
from renamers.renameWithReplicantCounts import renameWithReplicantCounts
from renamers.renameWithoutCount import removeCountsFromFastFile, removeCountsFromNamesFile
from renamers.updateNames import updateNames
from utils.joinFiles import joinFiles
from utils.splitKperFasta import splitK


def assemble_pear(args, debug=False):
    """Assembles reads from two (left and right) fastq files/directories.
    :param args: An argparse object with the following parameters:
                    name            Textual ID for the data set.
                    input_f         Forward Fastq Reads file or directory.
                    input_r         Reverse Fastq Reads file or directory.
                    threads         The number of threads to use durring assembly.
                    outdir          Directory where outputs will be saved.
    
    """
    # "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s -m %d"
    makeDirOrdie(args.outdir)
    forwards_reads = getInputFiles(args.input_f, "*_forward*", critical=False)
    forwards_reads += getInputFiles(args.input_f, "*_R1*", critical=False)
    reverse_reads = getInputFiles(args.input_r, "*_reverse*", critical=False)
    reverse_reads += getInputFiles(args.input_r, "*_R2*", critical=False)

    # Ensure that we have matching left and right reads
    if len(forwards_reads) != len(reverse_reads):
        print "Error: Unequal number of forwards/reverse reads."
        exit()

    if len(forwards_reads) == 0:
        print "Forwards reads should include the filename suffix \"_forward\" or \"R1\".  Reverse reads should \
                        include the filename suffix \"_reverse\" or \"R2\"."

    inputs = zip(set(forwards_reads), set(reverse_reads))
    pool = init_pool(min(len(inputs), args.threads))
    printVerbose("\tAssembling reads with pear")
    debugPrintInputInfo(inputs, "assemble")
    rslt = parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.ASSEMBLE_PEAR, [forwards, reverse,
                                    "%s/%s_%s" % (args.outdir, args.name, getFileName(forwards)), args.pearthreads],
                                    {"exists": [forwards, reverse], "positive": [args.pearthreads]})
                                    for forwards, reverse in inputs], pool)

    printVerbose("Done assembling sequences...")
    # Grab all the auxillary files (everything not containing ".assembled."
    aux_files = getInputFiles(args.outdir, "*", "*.assembled.*", ignore_empty_files=False)
    # make aux dir for extraneous files and move them there
    bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    cleanup_pool(pool)


def splitOnBarcodes(args, debug=False):
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
    makeDirOrdie(args.outdir)

    files_to_split = getInputFiles(args.input)

    file_id = range(len(files_to_split))
    file_id_pairs = zip(files_to_split, file_id)
    pool = init_pool(min(len(file_id_pairs), args.threads))
    printVerbose("Demuxing sequences...")
    debugPrintInputInfo(files_to_split, "demux")
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.DEMUX_FASTX,
                                                      [input_, args.barcodes, "%s/" % args.outdir,
                                                       "_%d_splitOut.fastq" % id_], {"exists": [input_]})
                                        for input_, id_ in file_id_pairs], pool)
    printVerbose("Demuxed sequences.")

    # gather output files and move them to their final destination
    output_files = enumerateDir(".", "*_splitOut_")
    bulk_move_to_dir(output_files, args.outdir)

    # Grab all the auxillary files
    aux_files = getInputFiles(args.outdir, "unmatched_*")
    # make aux dir for extraneous files and move them there
    bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    cleanup_pool(pool)


def renameSequences(args, debug=False):
    """Renames sequences in a fasta/fastq file as <filename>_ID0, <filename>_ID1, <filename>_ID2, etc., where <filename>
        is the name of the fasta/fastq file without any extensions or chewbacca suffixes.
    :param args: An argparse object with the following parameters:
                    input       Input file or folder containing only fasta or fastq files.
                    outdir      Directory where outputs will be saved.
                    filetype    File type of the input files.  Either 'fasta' or 'fastq'.
    """
    # Make the output directory, or abort if it already exists
    makeDirOrdie(args.outdir)

    # Gather input files
    inputs = getInputFiles(args.input)
    debugPrintInputInfo(inputs, "rename")
    pool = init_pool(min(len(inputs), args.threads))
    printVerbose("Renaming sequences...")
    # Run serialRename in parallel
    parallel(runPythonInstance,
             [(serialRename,
               input_, "%s/%s_renamed%s" % (args.outdir, strip_ixes(input_), os.path.splitext(input_)[1]),
               args.filetype, args.clip) for input_ in inputs], pool)
    printVerbose("Done renaming sequences...")

    # Grab all the auxillary files
    aux_files = getInputFiles(args.outdir, "*.names")
    # make aux dir for extraneous files and move them there
    bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    cleanup_pool(pool)


def trim_flexbar(args, debug=False):
    """Use Flexbar to trim adapters (and preceeding barcodes) from sequences in input file(s).  Sequences should be in
        the following format: <BARCODE><ADAPTER><SEQUENCE><RC_ADAPTER>, where ADAPTER is defined in the adapters file,
        and RC_ADAPTER is defined in the rcadapters file.  By default, Flexbar does not allow any 'N' characters in
        SEQUENCE, and will toss any sequences that do contain 'N'.  To avoid this, use the -u or --allowedns flags to
        specify the maximum number of 'N's to allow.
    # TODO fill this in.
    :param args:
    :param debug:
    """
    # "flexbar":  "flexbar -r \"%s\" -t \"%s\" -ae \"%s\" -a \"%s\"",
    makeDirOrdie(args.outdir)
    temp_file_name_template = "%s/temp_%s"
    debarcoded_file_name_template = "%s/%s_debarcoded"
    # TODO those exists validators dont really need ot be there since we globbed our files
    printVerbose("Trimming barcodes and adapters with flexbar")
    inputs = getInputFiles(args.input)
    pool = init_pool(min(len(inputs), args.threads))
    debugPrintInputInfo(inputs, "trim adapters from")
    # Trim the left
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.TRIM_FLEXBAR,
                                                      [input_file,
                                                       temp_file_name_template % ( args.outdir, strip_ixes(input_file)),
                                                       "LEFT", args.adapters, args.allowedns],
                                                      {"exists": [input_file, args.adapters]})
                                                        for input_file in  inputs], pool)

    temp_files = getInputFiles(args.outdir, "temp_*")
    debugPrintInputInfo(temp_files, "trim adapters from")

    # Trim the right
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.TRIM_FLEXBAR,
                                                      [input_file,
                                                       debarcoded_file_name_template % ( args.outdir,
                                                       strip_ixes(input_file)[5:]), "RIGHT", args.adaptersrc,
                                                       args.allowedns],
                                                      {"exists": [input_file, args.adaptersrc]})
                                                        for input_file in temp_files], pool)
    printVerbose("Done Trimming sequences.")

    # Grab all the auxillary files (everything not containing ".assembled."
    for file_ in getInputFiles(args.outdir, "temp_*"):
        os.remove(file_)

    cleanup_pool(pool)


def trimmomatic(args, debug=False):
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
    inputs = getInputFiles(args.input)
    pool = init_pool(min(len(inputs), args.threads))

    debugPrintInputInfo(inputs, "clean")
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.CLEAN_TRIMMOMATIC, [input_,
                                                    "%s/%s_cleaned.fastq" % (args.outdir, strip_ixes(input_)),
                                                    args.windowSize, args.quality, args.minLen],
                                                    {"exists": [args.outdir, input_],
                                                    "positive": [args.windowSize, args.quality, args.minLen]})
                                                    for input_ in inputs], pool)
    printVerbose("Done cleaning sequences.")

    cleanup_pool(pool)


def dereplicate(args, debug=False):
    makeDirOrdie(args.outdir)
    inputs = getInputFiles(args.input)
    pool = init_pool(min(len(inputs), args.threads))
    # REMOVES COUNTS FROM SEQUENCE NAMES IN ORDER TO CLUSTER PROPERLY
    # strip counts if we need to.
    stripCounts = args.stripcounts
    print stripCounts
    if stripCounts:
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

    # LOG DEREPLICATED SEQUENCES INTO A .NAMES FILE
    # generates a .names file named _uc_parsed.out
    # python parseUCtoNames.py uc.out uc_parsed.out
    input_ucs = getInputFiles(args.outdir, "*_uc.out")
    printVerbose("Generating a names file from dereplication.")
    debugPrintInputInfo(inputs, "parsed (into a names file)")
    parallel(runPythonInstance,
             [(parseUCtoNames, input_, "%s/%s_derep.names" % (args.outdir, strip_ixes(input_)))
              for input_ in input_ucs], pool)

    most_recent_names_files = getInputFiles(args.outdir, "*_derep.names")

    # UPDATE THE MOST CURRENT NAMES FILES WITH DEREPLICATION COUNTS
    if args.namesfile is not None:
        # Grab the old names file and the dereplicated names file
        old_names_files = getInputFiles(args.namesfile)
        derep_names_files = getInputFiles(args.outdir, "*_derep.names")

        printVerbose("Updating .names files with dereplicated data")
        logging.debug("%d Reference (old) names files to be read:" % len(old_names_files))
        logging.debug(str(old_names_files))
        logging.debug("%d Dereplicated (new) names files to be read:" % len(derep_names_files))
        logging.debug(str(derep_names_files))
        # updateNames (old_names_files, new_names_files, updated)
        updateNames(old_names_files, derep_names_files, args.outdir, "dereplicated")
        most_recent_names_files = getInputFiles(args.outdir, "dereplicated*")
        printVerbose("Done updating .names files.")

    if len(inputs) != len(most_recent_names_files):
        print ("Error: Number of input fastas (%d) is not equal to the number of names files (%d)." %
               (len(inputs), len(most_recent_names_files)))
        exit()
    fasta_names_pairs = zip(inputs, most_recent_names_files)
    # ADD COUNT TO SEQUENCE NAMES AND SORT BY COUNT
    # python renameWithReplicantCounts.py
    #               8_macse_out/MACSEOUT_MERGED.fasta uc_parsed.out dereplicated_renamed.fasta
    printVerbose("Adding dereplication data to unique fasta")
    parallel(runPythonInstance,
             [(renameWithReplicantCounts, fasta, names,
               "%s/%s_counts.fasta" % (args.outdir, strip_ixes(fasta)), 'fasta')
              for fasta, names in fasta_names_pairs], pool)
    printVerbose("Done adding data")

    aux_dir = makeAuxDir(args.outdir)
    names_dir = makeDirOrdie("%s_names_files" % args.outdir)
    bulk_move_to_dir(most_recent_names_files, names_dir)
    aux_files = getInputFiles(args.outdir, '*', "*_counts.fasta")
    bulk_move_to_dir(aux_files, aux_dir)


def partition(args, debug=False):
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
    inputs = getInputFiles(args.input)
    debugPrintInputInfo(inputs, "partitioned")
    printVerbose("Partitioning Files...")

    pool = init_pool(min(len(inputs), args.threads))

    parallel(runPythonInstance, [
        (splitK, input_, "%s/%s" % (args.outdir, strip_ixes(getFileName(input_))), args.chunksize, args.filetype)
        for input_ in inputs], pool)
    printVerbose("Done partitioning files.")
    cleanup_pool(pool)


def merge(args, debug=False):
    """Blindly concatenates files in a directory.
       :param args: An argparse object with the following parameters:
                       input        Cleaned inputs File.
                       outdir       Directory where outputs will be saved.
                       name         Name prefix for the merged file.
                       fileext      Output file extension.  e.g 'fasta', 'fastq', 'txt'
       """
    makeDirOrdie(args.outdir)
    inputs = getInputFiles(args.input)
    debugPrintInputInfo(inputs, "merged together.")
    printVerbose("Merging files...")
    joinFiles(inputs, "%s/%s_MERGED.%s" % (args.outdir, args.name, args.fileext))
    pool = init_pool(min(len(inputs), args.threads))
    printVerbose("Done merging.")
    cleanup_pool(pool)


def ungapFasta(args, debug=False):
    """Removes a gap character from sequences in a fasta file.  Useful for removing characters from an alignment file.
       :param args: An argparse object with the following parameters:
                       input        Fasta files to ungap.
                       outdir       Directory where outputs will be saved.
                       gapchar      A gap character to remove from the sequences.
                       fileext      Either 'fastq' or 'fasta'.
       """
    makeDirOrdie(args.outdir)
    input_files = getInputFiles(args.input)
    debugPrintInputInfo(input_files, "ungap.")
    pool = init_pool(min(len(input_files), args.threads))
    printVerbose("Removing all '%s' from sequences..." % args.gapchar)
    # ungap(file_to_clean, output_file_name, gap_char, file_type):
    parallel(runPythonInstance, [(ungap, input_, "%s/%s_cleaned.%s" % (args.outdir, strip_ixes(input_), 'fasta'),
                                  args.gapchar, args.fileext) for input_ in input_files], pool)
    printVerbose("Done removing.")
    cleanup_pool(pool)


def cluster(args, debug=False):
    # TODO what is the correct default behavior if no names file is supplied? Currently singletons.
    """Clusters sequences.
    :param args: An argparse object with the following parameters:
                    input       A file or folder containing fasta files to cluster.
                    output      The output directory results will be written to.
                    namesFile	A names file or folder containing names files that describe the input.
                                Note: if no names file is supplied, then entries in the fasta file are assumed to be
                                    singleton sequences. See <http://www.mothur.org/wiki/Name_file>
                    stripcounts If False, leaves any abbundance suffixes (_###) intact.  This may hurt clustering.
                                    Defaults to True.
    """

    # TODO IMPORTANT: make sure that the abundances in dereplicated_renamed.fasta are sorted in decreasing order
    # We need to explore this more. Run by default for now....
    # Nore that any program can be used here as long as two files
    # 1- seeds: contains the cluster centroids. This file contains the updates counts for each cluster.
    # ex. a seq 97_2 from the cluster, if selected as a seed, would not be for example 97_100. This indicates
    # that 98 sequences are now assigned to the cluster for which 97 is a seed.
    # 2-clustering.out: contains the clustering results. (see file for sample format)

    # ~/bin/swarm/src/swarm dereplicated_renamed.fasta \
    #  		      				       --output-file clustering.out -u uclust_file -w seeds
    # "swarm": program_paths["SWARM"] + " \"%s\" --output-file \"%s\" \
    #                                            -u \"%s\" -w \"%s\"",
    # CLUSTER

    makeDirOrdie(args.outdir)
    # Grab the fasta file(s) to cluster
    inputs = getInputFiles(args.input)
    pool = init_pool(min(len(inputs), args.threads))
    # Try to grab a names file
    most_recent_names_files = getInputFiles(args.namesfile, critical=False)
    debugPrintInputInfo(inputs, "clustered")
    if len(most_recent_names_files) != 0:
        debugPrintInputInfo(most_recent_names_files, "used as group listings")
        printVerbose("Input names files: %s\n" % most_recent_names_files)
    else:
        printVerbose("No name files provided, assuming singletons...\n")

    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.CLUSTER_SWARM, [input_,
                                                                "%s/%s_clustered.names" % (
                                                                    args.outdir, strip_ixes(input_)),
                                                                "%s/%s_clustered_uclust" % (
                                                                    args.outdir, strip_ixes(input_)),
                                                                "%s/%s_clustered_seeds" % (
                                                                    args.outdir, strip_ixes(input_))],
                                                      {"exists": [input_]}) for input_ in inputs], pool)

    # REMOVE COUNTS FROM CLUSTERING NAMES FILE
    # Grab the current names file and the new clustered names file (which needs to be cleaned)
    clustered_names_files = getInputFiles(args.outdir, "*_clustered.names")
    # Remove counts from the clustering names files
    printVerbose("Cleaning the .names file from clustering")
    debugPrintInputInfo(clustered_names_files, "cleaned")
    parallel(runPythonInstance,
             [(removeCountsFromNamesFile, input_, "%s/%s_uncount.names" % (args.outdir, strip_ixes(input_)))
              for input_ in clustered_names_files], pool)

    # UPDATE THE NAMES FILES WITH NEW CLUSTERS
    cleaned_clustered_names_files = getInputFiles(args.outdir, "*clustered_uncount.names")
    printVerbose("Updating .names files with clustering data")
    # updateNames (old_names_files, new_names_files, updated)
    updateNames(most_recent_names_files, cleaned_clustered_names_files, args.outdir, "postcluster")
    printVerbose("Done updating .names files.")

    printVerbose("Capitalizing sequences")
    # Convert the seeds files to uppercase (swarm writes in lowercase)
    inputs = getInputFiles(args.outdir, "*_seeds")
    parallel(runPythonInstance,
             [(capitalizeSeqs, input_, "%s.fasta" % input_) for input_ in inputs], pool)
    printVerbose("Done capitalizing sequences")

    printVerbose("Converting seeds files to .names files")
    # delete seeds file
    for input_ in inputs:
        os.remove(input_)
        inputs = getInputFiles(args.outdir, "*_seeds.fasta")
        parallel(runPythonInstance,
                 [(seedToNames, input_, "%s/%s.names" % (args.outdir, getFileName(input_))) for input_ in inputs],
                 pool)
    printVerbose("Done converting seeds files.")

    printVerbose("Moving aux files")
    # Gather and move the post-clustering names file
    names_dir = makeDirOrdie(args.outdir + "_names_files")
    post_cluster_names_file = getInputFiles(args.outdir, "postcluster_updated.names", critical=False)
    bulk_move_to_dir(post_cluster_names_file, names_dir)

    # Gather and move auxillary files
    aux_files = getInputFiles(args.outdir, "*", "*_seeds.fasta")
    bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    cleanup_pool(pool)


def query_fasta(args, debug=False):
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
    query_fastas = getInputFiles(args.input)
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
                                            [args.threads, query_fasta, ref_fasta,
                                            "%s/%s.out" % (args.outdir, strip_ixes(query_fasta)),
                                            "%s/%s.alnout" % (args.outdir, strip_ixes(query_fasta)), aln_user_string],
                                            {"exists": [query_fasta, ref_fasta], "positive":[args.threads]})
                                            for query_fasta, (ref_fasta, tax_info) in inputs], pool)
    printVerbose("Done querying.")
    printVerbose("Parsing output...")


    # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
    # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
    #
    # parseVSearchOutputAgainstFasta(vsearch_outfile, taxInfo, output_file, min_simmilarity, min_coverage):
    parallel(runPythonInstance, [(parseVSearchOutputAgainstFasta, "%s/%s.out" % (args.outdir, strip_ixes(query_fasta)),
                                            tax_info, "%s/%s_result.out" % (args.outdir, strip_ixes(query_fasta)),
                                            args.simmilarity, args.coverage)
                                            for query_fasta, (ref_fasta, tax_info) in inputs], pool)
    printVerbose("\nDone parsing...")
    # Gather and move auxillary files
    aux_files = getInputFiles(args.outdir, "*", "*_result.out")
    bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    cleanup_pool(pool)


def queryNCBI(args, debug=False):
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
    #		--userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt

    # expecting a fasta to annotate
    query_fastas = getInputFiles(args.input)
    debugPrintInputInfo(query_fastas, "queried agains NCBI.")
    pool = init_pool(min(len(query_fastas), args.threads))
    printVerbose("Querying NCBI...")
    # input, db, userout, alnout, aln_userfield
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.ALIGN_VSEARCH,
                                                        [args.threads, query_fasta, coi_fasta,
                                                        "%s/%s.out" % (args.outdir, strip_ixes(query_fasta)),
                                                        "%s/%s.alnout" % (args.outdir, strip_ixes(query_fasta)),
                                                        aln_user_string],
                                                        {"exists": [query_fasta, coi_fasta], "positive": [args.threads]})
                                                        for query_fasta in query_fastas], pool)

    printVerbose("Done with query.")
    printVerbose("Parsing output...")
    # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
    # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
    #
    # parseVSearchOutputAgainstNCBI(vsearch_out, ncbi_db, min_coverage, min_similarity)> parsed_nt.out
    parallel(runPythonInstance,
             [(parseVSearchOutputAgainstNCBI, "%s/%s.out" % (args.outdir, strip_ixes(query_fasta)), ncbi_db_string,
               "%s/%s_result.out" % (args.outdir, strip_ixes(query_fasta)),
               97, 85) for query_fasta in query_fastas], pool)
    printVerbose("Done processing.")

    # Gather and move auxillary files
    aux_dir = makeAuxDir(args.outdir)
    aux_files = getInputFiles(args.outdir, "*", "*_result.out")
    bulk_move_to_dir(aux_files, aux_dir)

    cleanup_pool(pool)


# TODO doc and debug
def build_matrix(args, debug=False):
    """Builds the unannotated OTU table.
    :param args: An argparse object with the following parameters:
                    names       File/folder with .names files.
                    outdir      Directory to put the output files.
                    groups      File/folder with .groups files.
                    barcodes    A single barcodes file listing all possible sample names.
    """
    makeDirOrdie(args.outdir)
    names = getInputFiles(args.names)
    debugPrintInputInfo(names, "read.")
    groups = getInputFiles(args.groups)
    debugPrintInputInfo(groups, "read.")
    barcodes = getInputFiles(args.barcodes)
    debugPrintInputInfo(barcodes, "read.")
    printVerbose("Building matrix...")
    buildOTUtable(names, groups, barcodes[0], "%s/%s.txt" % (args.outdir, "matrix"))
    printVerbose("Done building.")


# TODO doc and debug
def annotate_matrix(args, debug=False):
    """Annotates an OTU table.
    :param args: An argparse object with the following parameters:
                    input       File/folder with matrix files.
                    outdir      Directory to put the output files.
                    map      File/folder with .groups files.
    """
    makeDirOrdie(args.outdir)
    matricies = getInputFiles(args.input)
    debugPrintInputInfo(matricies, "annotated.")
    annotations = getInputFiles(args.annotation)
    debugPrintInputInfo(annotations, "parse.")
    inputs = product(matricies, annotations)

    pool = init_pool(min(len(matricies)*len(annotations), args.threads))
    printVerbose("Annotating matrix...")
    parallel(runPythonInstance, [(annotateOTUtable, matrix, annotation, "%s/%s.txt" % (args.outdir, "matrix"))
                                 for matrix, annotation in inputs], pool)
    printVerbose("Done Annotating.")

    cleanup_pool(pool)


def makeFasta(args, debug=False):
    """Converts a fastq file to fasta format.
    :param args: An argparse object with the following parameters:
                    input       File/folder with fastq files.
                    outdir      Directory to put the output files.xs
    """
    makeDirOrdie(args.outdir)
    inputs = getInputFiles(args.input)
    debugPrintInputInfo(inputs, "convert to fasta.")
    printVerbose("Converting to fasta...")

    pool = init_pool(min(len(inputs), args.threads))
    parallel(runPythonInstance, [(translateFastqToFasta, input_, "%s/%s.fasta" % (args.outdir, getFileName(input_)))
                                 for input_ in inputs], pool)
    printVerbose("Done converting.")

    cleanup_pool(pool)
