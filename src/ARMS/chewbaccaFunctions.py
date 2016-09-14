from itertools import product
from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from converters.capitalizeSeqs import capitalizeSeqs
from converters.fastqToFasta import translateFastqToFasta
from converters.seedToGroups import seedToGroups
from converters.ungap import ungap
from otu_tables.annotateOTUtable import annotateOTUtable
from otu_tables.buildOTUtable import buildOTUtable
from parsers.parseUCtoGroups import parseUCtoGroups
from parsers.parseVSearchoutForTaxa import parseVSearchOutputAgainstFasta, parseVSearchOutputAgainstNCBI
from renamers.renameSerially import serialRename
from renamers.renameWithReplicantCounts import renameWithReplicantCounts
from renamers.renameWithoutCount import removeCountsFromFastFile, removeCountsFromGroupsFile
from renamers.updateGroups import update_groups
from utils.joinFiles import joinFiles
from utils.splitKperFasta import splitK


def assemble_pear(args):
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
    # "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s -m %d"
    makeDirOrdie(args.outdir)
    forwards_reads = getInputFiles(args.input_f, "*_forward*", critical=False)
    reverse_reads = getInputFiles(args.input_r, "*_reverse*", critical=False)

    if len(forwards_reads) == 0 and len(reverse_reads) == 0:
        forwards_reads = getInputFiles(args.input_f, "*_R1*", critical=False)
        reverse_reads = getInputFiles(args.input_r, "*_R2*", critical=False)

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
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.ASSEMBLE_PEAR, [forwards, reverse,
                                                                                            "%s/%s_%s" % (
                                                                                                args.outdir, args.name,
                                                                                                getFileName(forwards)),
                                                                                            args.pearthreads],
                                                      {"exists": [forwards, reverse], "positive": [args.pearthreads]})
                                        for forwards, reverse in inputs], pool)

    printVerbose("Done assembling sequences...")
    # Grab all the auxillary files (everything not containing ".assembled."
    aux_files = getInputFiles(args.outdir, "*", "*.assembled.*", ignore_empty_files=False)
    # make aux dir for extraneous files and move them there
    bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    cleanup_pool(pool)


def split_on_barcodes(args):
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

    samples_dir = makeDirOrdie("%s_samples" % args.outdir)
    samples_files = getInputFiles(args.outdir, "*.samples")
    bulk_move_to_dir(samples_files, samples_dir)

    cleanup_pool(pool)


def trim_flexbar(args):
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
    # "flexbar":  "flexbar -r \"%s\" -t \"%s\" -ae \"%s\" -a \"%s\"",
    makeDirOrdie(args.outdir)
    temp_file_name_template = "%s/temp_%s"
    debarcoded_file_name_template = "%s/%s_debarcoded"
    # TODO those exists validators dont really need ot be there since we globbed our files
    printVerbose("Trimming barcodes and adapters with flexbar")
    inputs = getInputFiles(args.input)
    pool = init_pool(min(len(inputs), args.threads))
    debugPrintInputInfo(inputs, "trim adapters from")

    # Trim adapters from the left
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.TRIM_FLEXBAR,
                                                      [input_file,
                                                       temp_file_name_template % (args.outdir, strip_ixes(input_file)),
                                                       "LEFT", args.adapters, args.allowedns],
                                                      {"exists": [input_file, args.adapters]})
                                        for input_file in inputs], pool)

    temp_files = getInputFiles(args.outdir, "temp_*")
    debugPrintInputInfo(temp_files, "trim adapters from")

    # Trim the reverse complemented adapters from the right
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.TRIM_FLEXBAR,
                                                      [input_file,
                                                       debarcoded_file_name_template % (args.outdir,
                                                                                        strip_ixes(input_file)[5:]),
                                                       "RIGHT", args.adaptersrc,
                                                       args.allowedns],
                                                      {"exists": [input_file, args.adaptersrc]})
                                        for input_file in temp_files], pool)
    printVerbose("Done Trimming sequences.")

    # Grab all the auxillary files (everything not containing ".assembled."
    for file_ in getInputFiles(args.outdir, "temp_*"):
        os.remove(file_)

    cleanup_pool(pool)


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
    inputs = getInputFiles(args.input)
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
    inputs = getInputFiles(args.input)
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
    inputs = getInputFiles(args.input)
    debugPrintInputInfo(inputs, "partitioned")
    printVerbose("Partitioning Files...")

    pool = init_pool(min(len(inputs), args.threads))

    parallel(runPythonInstance, [
        (splitK, input_, "%s/%s" % (args.outdir, strip_ixes(getFileName(input_))), args.chunksize, args.filetype)
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
    inputs = getInputFiles(args.input)
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
    input_files = getInputFiles(args.input)
    debugPrintInputInfo(input_files, "ungap.")
    pool = init_pool(min(len(input_files), args.threads))
    printVerbose("Removing all '%s' from sequences..." % args.gapchar)
    # ungap(file_to_clean, output_file_name, gap_char, file_type):
    parallel(runPythonInstance, [(ungap, input_, "%s/%s_cleaned.%s" % (args.outdir, strip_ixes(input_), 'fasta'),
                                  args.gapchar, args.fileext) for input_ in input_files], pool)
    printVerbose("Done removing.")
    cleanup_pool(pool)


def cluster(args):
    """Clusters sequences.
    :param args: An argparse object with the following parameters:
                    input       A file or folder containing fasta files to cluster.
                    output      The output directory results will be written to.
                    gropusfile	Agroups file or folder containinggroups files that describe the input.
                                Note: if no groups file is supplied, then entries in the fasta file are assumed to be
                                    singleton sequences. See <http://www.mothur.org/wiki/Name_file>
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
    # Try to grab agroups file
    most_recent_groups_files = getInputFiles(args.groupsfile, critical=False)
    debugPrintInputInfo(inputs, "clustered")
    if len(most_recent_groups_files) != 0:
        debugPrintInputInfo(most_recent_groups_files, "used as sample listings")
        printVerbose("Inputgroups files: %s\n" % most_recent_groups_files)
    else:
        printVerbose("No name files provided, assuming singletons...\n")

    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.CLUSTER_SWARM, [input_,
                                                                                            "%s/%s_clustered.groups" % (
                                                                                                args.outdir,
                                                                                                strip_ixes(input_)),
                                                                                            "%s/%s_clustered_uclust" % (
                                                                                                args.outdir,
                                                                                                strip_ixes(input_)),
                                                                                            "%s/%s_clustered_seeds" % (
                                                                                                args.outdir,
                                                                                                strip_ixes(input_))],
                                                      {"exists": [input_]}) for input_ in inputs], pool)

    # REMOVE COUNTS FROM CLUSTERING GROUPS FILE
    # Grab the current groups file and the new clustered groups file (which needs to be cleaned)
    clustered_groups_files = getInputFiles(args.outdir, "*_clustered.groups")
    # Remove counts from the clustering groups files
    printVerbose("Cleaning the .groups file from clustering")
    debugPrintInputInfo(clustered_groups_files, "cleaned")
    parallel(runPythonInstance,
             [(removeCountsFromGroupsFile, input_, "%s/%s_uncount.groups" % (args.outdir, strip_ixes(input_)))
              for input_ in clustered_groups_files], pool)

    # UPDATE THE NAMES FILES WITH NEW CLUSTERS
    cleaned_clustered_groups_files = getInputFiles(args.outdir, "*clustered_uncount.groups")
    printVerbose("Updating .groups files with clustering data")
    # update_groups (old_groups_files, new_groups_files, updated)
    update_groups(most_recent_groups_files, cleaned_clustered_groups_files, args.outdir, "postcluster")
    printVerbose("Done updating .groups files.")

    printVerbose("Capitalizing sequences")
    # Convert the seeds files to uppercase (swarm writes in lowercase)
    inputs = getInputFiles(args.outdir, "*_seeds")
    parallel(runPythonInstance,
             [(capitalizeSeqs, input_, "%s.fasta" % input_) for input_ in inputs], pool)
    printVerbose("Done capitalizing sequences")

    printVerbose("Converting seeds files to .groups files")
    # delete seeds file
    for input_ in inputs:
        os.remove(input_)
        inputs = getInputFiles(args.outdir, "*_seeds.fasta")
        parallel(runPythonInstance,
                 [(seedToGroups, input_, "%s/%s.groups" % (args.outdir, getFileName(input_))) for input_ in inputs],
                 pool)
    printVerbose("Done converting seeds files.")

    printVerbose("Moving aux files")
    # Gather and move the post-clusteringgroups file
    groups_dir = makeDirOrdie(args.outdir + "_groups_files")
    post_cluster_groups_file = getInputFiles(args.outdir, "postcluster_updated.groups", critical=False)
    bulk_move_to_dir(post_cluster_groups_file, groups_dir)

    # Gather and move auxillary files
    aux_files = getInputFiles(args.outdir, "*", "*_seeds.fasta")
    bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    cleanup_pool(pool)


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
    query_fastas = getInputFiles(args.input)
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
    samples = getInputFiles(args.groups)
    debugPrintInputInfo(samples, "read.")
    barcodes = getInputFiles(args.barcodes)
    debugPrintInputInfo(barcodes, "read.")
    printVerbose("Building matrix...")
    buildOTUtable(groups, samples, barcodes[0], "%s/%s.txt" % (args.outdir, "matrix"))
    printVerbose("Done building.")


# TODO doc and debug
def annotate_matrix(args):
    """Annotates an OTU table.
    :param args: An argparse object with the following parameters:
                    input       File/folder with matrix files.
                    outdir      Directory to put the output files.
                    map      File/folder with .samples files.
    """
    makeDirOrdie(args.outdir)
    matricies = getInputFiles(args.input)
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
    inputs = getInputFiles(args.input)
    debugPrintInputInfo(inputs, "convert to fasta.")
    printVerbose("Converting to fasta...")

    pool = init_pool(min(len(inputs), args.threads))
    parallel(runPythonInstance, [(translateFastqToFasta, input_, "%s/%s.fasta" % (args.outdir, getFileName(input_)))
                                 for input_ in inputs], pool)
    printVerbose("Done converting.")

    cleanup_pool(pool)
