import subprocess
from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner
from converters.capitalizeSeqs import capitalizeSeqs
from converters.fastqToFasta import translateFastqToFasta
from converters.seedToNames import seedToNames
from filters.prescreen import screen
from filters.removeMacseRefs import removeMacseRefs
from filters.ungap import ungap
from parsers.getSeedSequences import getSeedSequences
from parsers.getTaxFromId import parseVSearchToTaxa
from parsers.parseVSearchout import parseVSearchout
from renamers.formatWithSwarmCounts import formatWithSwarmCounts
from renamers.renameSequences import serialRename
from renamers.renameWithCount import renameSequencesWithCount
from utils.joinFiles import joinFiles
from utils.splitKperFasta import splitK


def assemble(args, pool=Pool(processes=1)):
    """Assembles reads from two (left and right) fastq files.  Handler for the mothur and pear assemblers.  Validates
        argument dependencies.
    :param args: An argparse object with the following parameters:
                program         The program to use.  Either 'pear' or 'mothur'
                name		    Assembled File Prefix.
                input_f		    Forward Fastq Reads file or directory.
                input_r		    Reverse Fastq Reads file or directory.
                outdir		    Directory where outputs will be saved.

                # options for mothur
                oligos          Oligos file with barcode and primer sequences.
                bdiffs          # of allowed barcode mismatches.
                pdiffs          # of allowed primer mismatches.
                # options for pear
                threads         The number of threads to use per process. [Hint: threads * subprocesses <= # of processors on your machine]
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        supported_programs = {
            "mothur": assemble_mothur,
            "pear": assemble_pear
        }
        keys = supported_programs.keys()
        if args.program in keys:
            # Validate argument dependencies
            if args.program == "mothur":
                if not args.oligos:
                    "-g parameter required for mothur assembler"
            # Make the output directory
            makeDir(args.outdir)
            # Run and get a list of output files
            outputs = supported_programs[args.program](args, pool)

        else:
            "Invalid program choice.  Supported programs are " + str(keys)
    except KeyboardInterrupt:
        cleanupPool(pool)


def assemble_pear(args, pool=Pool(processes=1)):
    """Assembles reads from two (left and right) fastq files/directories.
    :param args: An argparse object with the following parameters:
                    name            Textual ID for the data set.
                    input_f         Forward Fastq Reads file or directory.
                    input_r         Reverse Fastq Reads file or directory.
                    threads         The number of threads to use durring assembly.
                    outdir          Directory where outputs will be saved.
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """


    # "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s -m %d"
    try:

        name_error_text = "Forwards reads should include the filename suffix \"_forward\" " \
                          + "or \"R1\".  Reverse reads should include the filename suffix \"_reverse\" or \"R2\"."

        forwards_reads = getInputs(args.input_f, "*_forward*")
        forwards_reads += getInputs(args.input_f, "*_R1*")
        reverse_reads = getInputs(args.input_r, "*_reverse*")
        reverse_reads += getInputs(args.input_r, "*_R2*")

        # Ensure that we have matching left and right reads
        if len(forwards_reads) != len(reverse_reads):
            print "Error: Unequal number of forwards/reverse reads."
            exit()

        inputs = zip(set(forwards_reads), set(reverse_reads))

        logging.debug("%d Input files to assemble:" % len(inputs))
        logging.debug(str(inputs))

        printVerbose("\tAssembling reads with pear")
        parallel(runProgramRunner, [ProgramRunner("pear", [forwards, reverse,
                                                           "%s/%s_%s" % (args.outdir, args.name, getFileName(forwards)),
                                                           args.threads],
                                                  {"exists": [forwards, reverse]})
                                    for forwards, reverse in inputs], pool)

        printVerbose("Done assembling sequences...")

        # Grab all the auxillary files (everything not containing ".assembled."
        aux_files = getInputs(args.outdir, "*", "*.assembled.*", ignoreEmptyFiles=False)
        # make aux dir for extraneous files and move them there
        bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    except KeyboardInterrupt:
        cleanupPool(pool)


#TODO handle file cleanup
# TODO handle debug/verbose printing
def assemble_mothur(args, pool=Pool(processes=1)):
    """Assembles reads with Mothur's make.contigs command.
    :param args: An argparse object with the following parameters:
                forward	    Forward read fastq file or folder.
                reverse	    Reverse read fastq file or folder.
                oligos	    Oligos file with barcode and primer sequences.  See <http://www.mothur.org/wiki/Oligos_File>
                bdiffs	    # of allowed barcode mismatches
                pdiffs 	    # of allowed primer mismatches
                procs	    Number of processors to use
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # "make.contigs": "mothur \'#make.contigs(ffastq=%s, rfastq=%s, bdiffs=1, pdiffs=2, oligos=%s, processors=%s)\'"
    try:
        printVerbose("\tAssembling reads with pear")
        parallel(runProgramRunner, [ProgramRunner("trimmomatic", [args.forward, args.reverse, args.bdiffs, args.pdiffs,
                                                                  args.oligos, args.processors],
                                                  {"exists": [args.outputFile, args.inputFile, args.oligos],
                                                   "positive": [args.processors],
                                                   "non-Negative": [args.bdiffs, args.pdiffs]}, pool)
                                    ])
        # TODO: move output files to outdir
    except KeyboardInterrupt:
        cleanupPool(pool)


def splitOnBarcodes(args, pool=Pool(processes=1)):
    """Splits a fasta/fastq file on a set of barcodes.  For a set of k input files, each input file is assigned a file_id.
    Then, each file is split into a set of subfiles, where each subfile named <sample>_splitOut_<k> contains the sequences belonging to <sample> from file <k>.
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
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # TODO MAHDI  Comment the rest of the program to replace the pipeline with
    # 1- split_libraries_fastq.py
    # 2- usearch -fastq_filter
    # "barcode.splitter": "cat \"%s\" | " + programPaths["FASTX"] + 'fastx_barcode_splitter.pl  --bcfile "%s" \
    #                                    -prefix "%s" --suffix %s --bol --mismatches 1',
    try:
        makeDir(args.outdir)

        files_to_split = getInputs(args.input)
        logging.debug("%d Input files to demux:" % len(files_to_split))
        logging.debug(files_to_split)

        file_id = range(len(files_to_split))
        file_id_pairs = zip(files_to_split, file_id)
        printVerbose("Demuxing sequences...")
        parallel(runProgramRunner, [ProgramRunner("barcode.splitter",
                                                  [input_, args.barcodes, "%s/" % args.outdir,
                                                   "_%d_splitOut.fastq" % id_],
                                                  {"exists": [input_]}) for input_, id_ in file_id_pairs], pool)
        printVerbose("Demuxed sequences.")

        # gather output files and move them to their final destination
        output_files = enumerateDir(".", "*_splitOut_")
        bulk_move_to_dir(output_files, args.outdir)

        # Grab all the auxillary files
        aux_files = getInputs(args.outdir, "unmatched_*")
        # make aux dir for extraneous files and move them there
        bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    except KeyboardInterrupt:
        cleanupPool(pool)


def renameSequences(args, pool=Pool(processes=1)):
    """Renames sequences in a fasta/fastq file as <filename>_ID0, <filename>_ID1, <filename>_ID2, etc., where <filename>
        is the name of the fasta/fastq file without any extensions or chewbacca suffixes.
    :param args: An argparse object with the following parameters:
                    input       Input file or folder containing only fasta or fastq files.
                    outdir      Directory where outputs will be saved.
                    filetype    File type of the input files.  Either 'fasta' or 'fastq'.
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        # Make the output directory, or abort if it already exists
        makeDir(args.outdir)

        # Gather input files
        inputs = getInputs(args.input)
        logging.debug("%d Input files to rename:" % len(inputs))
        logging.debug(str(inputs))

        printVerbose("Renaming sequences...")
        # Run serialRename in parallel
        parallel(runPython,
                 [(serialRename,
                   input_, "%s/%s_renamed%s" % (args.outdir, strip_ixes(getFileName(input_)),
                                                os.path.splitext(input_)[1]), args.filetype) for input_ in inputs],
                 pool)
        printVerbose("Done renaming sequences...")

        # Grab all the auxillary files
        aux_files = getInputs(args.outdir, "*.names")
        # make aux dir for extraneous files and move them there
        bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    except KeyboardInterrupt:
        cleanupPool(pool)


def trim(args, pool=Pool(processes=1)):
    """Trims the adapter and barcode from each sequence in a given input fasta file, or directory of fasta files.

    :param args: An argparse object with the following parameters:
            input       File or directory containing fasta files to trim.
            outdir      Directory where outputs will be saved.
            program         The program to use.  Either 'pear' or 'mothur'.

            # options for mothur
            oligos          Oligos file with barcode and primer sequences.  See <http://www.mothur.org/wiki/Oligos_File>
            # options for flexbar
            barcodes    Filepath to a barcodes file (Two column, tab delimited file where the first column contains the
                            barcode sequence, and the second column contains barcode names).
            adapters    Filepath to an adapters file (twocolumn
    :param pool:
    :return:
    """
    try:
        supported_programs = {
            "flexbar": trim_flexbar,
            "mothur": trim_mothur}
        keys = supported_programs.keys()
        prog = args.program
        if args.program in keys:
            # Validate argument dependencies
            if prog == "flexbar":
                if not args.barcodes:
                    "-b parameter requried for flexbar trim"
                if not args.adapters:
                    "-a parameter required for flexbar trim"

            if prog == "mothur":
                if not args.oligos:
                    "-g parameter required for mothur trim"

            # Make the output directory
            makeDir(args.outdir)
            # Run and get a list of output files
            supported_programs[prog](args, pool)
        else:
            print "Invalid program choice.  Supported programs are " + str(keys)
    except KeyboardInterrupt:
        cleanupPool(pool)


def trim_flexbar(args, pool=Pool(processes=1)):
    # "flexbar":  "flexbar -r \"%s\" -t \"%s\" -ae \"%s\" -a \"%s\"",
    try:
        makeDir(args.outdir)

        input_files = getInputs(args.input)
        logging.debug("%d InputFiles to trim adapters from:" % len(input_files))
        logging.debug(str(input_files))

        temp_file_name_template = "%s/temp_%s"
        debarcoded_file_name_template = "%s/%s_debarcoded"


        # TODO those exists validators dont really need ot be there since we globbed our files
        # TODO these should write to the input directory, then return the file names, so the parent can move and rename
        printVerbose("Trimming barcodes and adapters with flexbar")
        # Trim the left
        parallel(runProgramRunner, [ProgramRunner("flexbar",
                                          [input_file, temp_file_name_template % (args.outdir, strip_ixes(input_file)),
                                           "LEFT", args.barcodes],
                                          {"exists": [input_file]}) for input_file in input_files], pool)

        temp_files = getInputs(args.outdir, "temp_*")
        logging.debug("%d InputFiles to debarcode:" % len(temp_files))
        logging.debug(str(temp_files))

        # Trim the right
        parallel(runProgramRunner, [ProgramRunner("flexbar",
                                          [input_file,
                                           debarcoded_file_name_template % (args.outdir, strip_ixes(input_file)[5:]),
                                           "RIGHT", args.adapters],
                                          {"exists": [input_file]}) for input_file in temp_files], pool)
        printVerbose("Done Trimming sequences.")

        # Grab all the auxillary files (everything not containing ".assembled."
        for file_ in getInputs(args.outdir, "temp_*"):
            os.remove(file_)

    except KeyboardInterrupt:
        cleanupPool(pool)

# TODO debug/verbose printing
# TODO file cleanup
def trim_mothur(args, pool=Pool(processes=1)):
    """
    :param args: An argparse object with the following parameters:
                    inputFasta  Fasta file with sequences to be trimmed
                    oligos      A mothur oligos file with barcodes and primers.  See:
                                    <http://www.mothur.org/wiki/Trim.seqs#oligos>
                    outdir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        printVerbose("Trimming barcodes and adapters with mothur")
        inputs = getInputs(args.input)
        parallel(runProgramRunner, [ProgramRunner("trim.seqs", [input_, args.oligos],
                                                  {"exists": [args.oligos]})
                                    for input_ in inputs], pool)
        printVerbose("Trimmed sequences.")
    except KeyboardInterrupt:
        cleanupPool(pool)


def trimmomatic(args, pool=Pool(processes=1)):
    """Uses a sliding window to identify and trim away areas of low quality.

    :param args: An argparse object with the following parameters:
                    inputFile	Input Fastq file
                    outputFile	Output Fastq file
                    windowSize	Width of the sliding window
                    quality 	Minimum passing quality for the sliding window
                    minLen	    Minimum passing length for a cleaned sequence
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # "trimomatic":       "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
    # -%phred %input %output SLIDINGWINDOW:%windowsize:%minAvgQuality MINLEN:%minLen"
    try:
        makeDir(args.outdir)

        inputs = getInputs(args.input)
        logging.debug("%d Input files to clean:" % len(inputs))
        logging.debug(str(inputs))

        printVerbose("Cleaning sequences with Trimmomatic...")
        parallel(runProgramRunner, [ProgramRunner("trimmomatic", [input_,
                                                                  "%s/%s_cleaned.fastq" % (
                                                                      args.outdir, strip_ixes(input_)),
                                                                  args.windowSize, args.quality, args.minLen],
                                                  {"exists": [args.outdir, input_],
                                                   "positive": [args.windowSize, args.quality, args.minLen]})
                                    for input_ in inputs], pool)
        printVerbose("Done cleaning sequences.")

    except KeyboardInterrupt:
        cleanupPool(pool)

# TODO We should probably move the fasta conversion out to a different step
def dereplicate(args, pool=Pool(processes=1)):
    """Uses a sliding window to identify and trim away areas of low quality.

       :param args: An argparse object with the following parameters:
                       inputFile	Input Fastq file
                       outdir	    Output Fastq file
       :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
       """
    # ls *cleaned.fasta  | parallel  "/ARMS/bin/usearch7.0.1090 -derep_fulllength {} -output {/.}_derep.fa
    #                                           -uc {/.}_uc.out"
    try:
        makeDir(args.outdir)
        #TODO debug/verbose for conversion
        # grab all the _cleaned files and turn the fastqs into fastas
        input_fqs = getInputs(args.input)
        output_fas = []
        fastq_exts = set([".fq", ".fastq"])

        # Generate a list of files to translate from fastq to fasta
        for input_ in input_fqs:
            if os.path.splitext(input_)[1].lower() in fastq_exts:
                input_dir = os.path.dirname(input_)
                base_name = os.path.basename(input_)
                file_name = os.path.splitext(base_name)[0]
                fasta_name = "%s/%s.fa" % (args.outdir, file_name)
                translateFastqToFasta(input_, fasta_name)
                output_fas.append(fasta_name)

        # zip those names to parallel array s
        file_args = zip(input_fqs, output_fas)
        # Translate fastq files in parallel
        parallel(runPython, [(translateFastqToFasta, input_fq, output_fa) for input_fq, output_fa in file_args], pool)


        logging.debug("%d Input files to search for replication:" % len(output_fas))
        logging.debug(str(output_fas))
        printVerbose("Searching for identical reads...")
        # parse the number of identical reads included in each sequence and write them to the
        #       {sample_file}_parsed.out
        # ls *cleaned.fasta  | parallel  "~/ARMS/bin/usearch7.0.1090 -derep_fulllength {} -output {/.}_derep.fa
        #                                                   -uc {/.}_uc.out"
        parallel(runProgramRunner, [ProgramRunner("usearch", [input_,
                                                              "%s/%s_derep.fa" % (
                                                                  args.outdir, strip_ixes(getFileName(input_))),
                                                              "%s/%s_uc.out" % (
                                                                  args.outdir, strip_ixes(getFileName(input_)))],
                                                  {"exists": [args.outdir, input_]})
                                    for input_ in output_fas], pool)
        printVerbose("Done searching.")

        # Collect .uc files
        input_ucs = getInputs(args.outdir, "*_uc.out")
        logging.debug("%d Input .out files to parse:" % len(input_ucs))
        logging.debug(str(input_ucs))
        printVerbose("Parsing search logs for seed sequences...")
        # Collect seeds
        # #ls *uc.out | parallel "python ~/ARMS/bin/getSeedSequences.py {} {.}_parsed.out"
        parallel(runPython,
                 [(getSeedSequences, input_, "%s/%s.names" % (args.outdir, strip_ixes(getFileName(input_)))) for
                  input_ in input_ucs], pool)
        printVerbose("Done parsing search logs.")

        # Rename the sequences to include the the number of identical reads.
        #  Ex. in the fasta file {sample_file}_derep_renamed.fa, read 123_10 indicate that for sequence
        # which id is 123, there 10 sequences that identical to it and which were discarded.
        # ls * cleaned.fasta | parallel "python ~/ARMS/bin/renameWithCount.py {/.}_derep.fa {/.}_uc_parsed.out {/.}_derep_renamed.fa"
        input_fa = getInputs(args.outdir, "*_derep.fa")
        input_uc_parsed = getInputs(args.outdir, "*.names")
        inputs = zip(input_fa, input_uc_parsed)
        logging.debug("%d Input fasta to rename and countfile:" % len(inputs))
        logging.debug(str(inputs))

        # renameSequencesWithCount(input_fasta, count_file, outfile):
        printVerbose("Renaming sequences with replication counts...")
        parallel(runPython,
                 [(renameSequencesWithCount,
                   fasta_file, count_file,
                   "%s/%s_derepCount.fa" % (args.outdir, strip_ixes(getFileName(fasta_file))))
                  for fasta_file, count_file in inputs], pool)
        printVerbose("Done renaming sequences.")

        aux_files = getInputs(args.outdir, "*", "*_derepCount.*")
        bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))
    except KeyboardInterrupt:
        cleanupPool(pool)



def align_mothur(args, pool=Pool(processes=1)):
    """Aligns sequences with mothur

    :param args: An argparse object with the following parameters:
                input        Database against which to align and filter reads
                ref          Directory containig the samples to be cleaned
                outdir       Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    :return:
    """
    try:
        # Make the output directory, or abort if it already exists

        makeDir(args.outdir)

        inputs = getInputs(args.input)
        logging.debug("%d Input files to be aligned:" % len(inputs))
        logging.debug(str(inputs))

        printVerbose("Aligning sequences with mothur against BIOCODE")
        # ls *debarcoded_cleaned_derep_renamed.fa | parallel  \
        #   mothur '"#align.seqs(candidate={}, template=~/ARMS/data/BIOCODETEMPLATE, flip=t)"'
        #
        # "align.seqs": "mothur \'#align.seqs(candidate=%s, template=%s, flip=t)\'",
        parallel(runProgramRunner, [ProgramRunner("align.seqs", [input_, args.ref],
                                                  {"exists": [args.outdir, args.ref, input_]})
                                    for input_ in inputs], pool)
        printVerbose("Done with alignment.")

        out_exts = ["*.align", "*.align.report", "*.accnos"]

        # collect and move output files
        if len(inputs) > 0:
            indir = os.path.dirname(inputs[0])
            out_files = []
            for ext in out_exts:
                out_files += getInputs(indir, ext)
            out_files += getInputs(".", "mothur.*.logfile")
            out_files += getInputs(".", "8mer")
            bulk_move_to_dir(out_files, args.outdir)

        # rename output files
        for file_ in getInputs(args.outdir, "*"):
            os.rename(file_, "%s/%s%s" % (args.outdir, strip_ixes(file_), os.path.splitext(file_)[1]))

        # Grab the auxillary files and move them to an aux directory
        aux_files = getInputs(args.outdir, "*", "*.align")
        bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    except KeyboardInterrupt:
        cleanupPool(pool)


def partition(args, pool=Pool(processes=1)):
    """ Partition the cleaned file into chunks of user-defined size.

    :param args: An argparse object with the following parameters:
                    input	    Input fasta file to split
                    outdir      Directory where outputs will be saved
                    chunksize	Chunksize.
                    filetype	Filetype of the files to be partitioned
    """
    # def splitK(inputFasta, prefix, nbSeqsPerFile, filetype):
    try:
        # Make the output directory, or abort if it already exists
        makeDir(args.outdir)

        # Gather input files
        inputs = getInputs(args.input)

        # Run renamer in parallel
        parallel(runPython, [
            (splitK, input_, "%s/%s" % (args.outdir, strip_ixes(getFileName(input_))), args.chunksize, args.filetype)
            for input_ in inputs], pool)


    except KeyboardInterrupt:
        cleanupPool(pool)


def merge(args, pool=Pool(processes=1)):
    """Clusters sequences.
       :param args: An argparse object with the following parameters:
                       input        Cleaned inputs File.
                       outdir       Directory where outputs will be saved.
                       name         Name prefix for the merged file.
                       fileext      Output file extension.  e.g 'fasta', 'fastq', 'txt'
       :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
       """
    try:
        makeDir(args.outdir)
        inputs = getInputs(args.input)
        joinFiles(inputs, "%s/%s_MERGED.%s" % (args.outdir, args.name, args.fileext))

    except KeyboardInterrupt:
        cleanupPool(pool)

def ungapFasta(args, pool=Pool(processes=1)):
    """Clusters sequences.
       :param args: An argparse object with the following parameters:
                       input        Fasta files to ungap.
                       outdir       Directory where outputs will be saved.
                       gapchar      A gap character to remove from the sequences.
       :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
       """
    try:
        makeDir(args.outdir)
        inputs = getInputs(args.input)

        #ungap(file_to_clean, output_file_name, gap_char, file_type):
        parallel(runPython, [(ungap, input_, "%s/%s_cleaned.%s" % (args.outdir, strip_ixes(input_), 'fasta'),
                              args.gapchar, args.fileext)
                             for input_ in inputs], pool)
    except KeyboardInterrupt:
        cleanupPool(pool)


def cluster(args, pool=Pool(processes=1)):
    """Clusters sequences.
    :param args: An argparse object with the following parameters:
                    inputFasta	Cleaned inputs File
                    program     Program for detecting and removing chimeras. Default is uchime
                    ...and exactly one of the following:
                    namesFile	Reference .names file. See <http://www.mothur.org/wiki/Name_file>
                    refDB       Reference database file.
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        makeDir(args.outdir)
        ## prepare the file for clustering. Dereplicates across smaples and renames the resulting sequences with hashes
        # "vsearch": program_paths["VSEARCH"] + "--derep_fulllength \"%s\" --sizeout --fasta_width 0  \
        #                                    --output amplicons_linearized_dereplicated.fasta -uc ",
        inputs = getInputs(args.input)

        parallel(runProgramRunner, [ProgramRunner("vsearch.derep",
                                                  [input_, "%s/%s.fa" % (args.outdir, getFileName(input_)),
                                                   "%s/%s_uc.out" % (args.outdir, getFileName(input_))],
                                                  {"exists": [input_]}) for input_ in inputs], pool)

        # python getSeedSequences.py uc.out uc_parsed.out
        input_ucs = getInputs(args.outdir, "*_uc.out")
        parallel(runPython,
                 [(getSeedSequences, input_, "%s/%s_parsed.out" % (args.outdir, getFileName(input_))) for
                  input_ in input_ucs], pool)

        #TODO update the names file with updateNames.py


        # python formatWithSwarmCounts.py
        #                                   8_macse_out/MACSEOUT_MERGED.fasta uc_parsed.out dereplicated_renamed.fasta
        parallel(runPython,
                 [(formatWithSwarmCounts, input_, "%s/%s_uc_parsed.out" % (args.outdir, getFileName(input_)),
                   "%s/%s_derep_renamed.fasta" % (args.outdir, getFileName(input_))) for
                  input_ in inputs], pool)
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
        inputs = getInputs(args.outdir, "*_derep_renamed.fasta")

        parallel(runProgramRunner, [ProgramRunner("swarm", [input_,
                                                            "%s/%s_clustering.out" % (args.outdir, strip_ixes(input_)),
                                                            "%s/%s_uclust" % (args.outdir, strip_ixes(input_)),
                                                            "%s/%s_seeds" % (args.outdir, strip_ixes(input_))],
                                                  {"exists": [input_]}) for input_ in inputs], pool)
        ## We convert the seeds files to uppercase (strangely, it is output in lowercase by swarm)
        ## this might not be necessary with other clustering programs

        # write upppercase seeds fasta
        inputs = getInputs(args.outdir, "*_seeds")
        parallel(runPython,
                 [(capitalizeSeqs, input_, "%s.fasta" % input_) for input_ in inputs], pool)
        # delete seeds file
        for input_ in inputs:
            os.remove(input_)
            inputs = getInputs(args.outdir, "*_seeds.fasta")
            parallel(runPython,
                     [(seedToNames, input_, "%s/%s.names" % (args.outdir, getFileName(input_))) for input_ in inputs],
                     pool)
            # generate names file

        # Gather and move auxillary files
        aux_files = getInputs(args.outdir, "*", "*_seeds.fasta")
        bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))
    except KeyboardInterrupt:
        cleanupPool(pool)


def queryBiocode(args, pool=Pool(processes=1)):
    """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

    :param args: An argparse object with the following parameters:
                    accnosFile  List of sequence names to remove
                    outdir      Directory to put the output files
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        DB_STRING = os.path.expanduser("~/ARMS/data/BiocodePASSED_SAP.txt")
        ALN_USER_STRING = ""
        makeDir(args.outdir)

        # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
        #		--userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt

        # expecting a fasta to annotate
        inputs = getInputs(args.input)
        # input, db, userout, alnout, aln_userfield
        parallel(runProgramRunner, [ProgramRunner("vsearch.usearch_global",
                                                  [input_, DB_STRING, "%s/%s.out" % (args.outdir, strip_ixes(input_)),
                                                   "%s/%s.alnout" % (args.outdir, strip_ixes(input_)), ALN_USER_STRING],
                                                  {"exists": [input_]}) for input_ in inputs], pool)
        # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
        # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
        #
        # parseVSearchout.py ~/ARMS/data/BiocodePASSED_SAP_tax_info.txt vsearch_out  parsed_BIOCODE.out 97 85
        parallel(runPython, [(parseVSearchout, os.path.expanduser("~/ARMS/data/BiocodePASSED_SAP_tax_info.txt"),
                              "%s/%s.out" % (args.outdir, strip_ixes(input_)),
                              "%s/%s_parsed_BIOCODE.out" % (args.outdir, strip_ixes(input_)),
                              97, 85) for input_ in inputs], pool)

    except KeyboardInterrupt:
        cleanupPool(pool)


def queryNCBI(args, pool=Pool(processes=1)):
    """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

    :param args: An argparse object with the following parameters:
                    accnosFile  List of sequence names to remove
                    outdir      Directory to put the output files
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        COI_DB_STRING = os.path.expanduser("~/ARMS/refs/COI.fasta")
        NCBI_DB_STRING = os.path.expanduser("~/ARMS/refs/ncbi.db")
        ALN_USER_STRING = "--userfields query+target+id+alnlen+qcov"
        makeDir(args.outdir)

        # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
        #		--userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt

        # expecting a fasta to annotate
        inputs = getInputs(args.input)
        # input, db, userout, alnout, aln_userfield
        parallel(runProgramRunner, [ProgramRunner("vsearch.usearch_global",
                                                  [input_, COI_DB_STRING,
                                                   "%s/%s.out" % (args.outdir, strip_ixes(input_)),
                                                   "%s/%s.alnout" % (args.outdir, strip_ixes(input_)), ALN_USER_STRING],
                                                  {"exists": [input_]}) for input_ in inputs], pool)
        # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
        # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
        #
        # parseVSearchToTaxa(vsearch_out, ncbi_db, min_coverage, min_similarity)> parsed_nt.out
        parallel(runPython, [(parseVSearchToTaxa, "%s/%s.out" % (args.outdir, strip_ixes(input_)), NCBI_DB_STRING,
                              "%s/%s_parsed_nt.out" % (args.outdir, strip_ixes(input_)),
                              97, 85) for input_ in inputs], pool)

    except KeyboardInterrupt:
        cleanupPool(pool)


def makeFastq(args, pool=Pool(processes=1)):
    """Finds chimeric sequences from a fasta file and writes them to an accons file.
    :param args: An argparse object with the following parameters:
                    inputFasta	Input Fasta file
                    inputQual	Input Qual file
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        parallel(runProgramRunner, [ProgramRunner("make.fastq", [args.inputFasta, args.inputQual],
                                                  {"exists": [args.inputFasta, args.inputQual]}, pool)
                                    ])

    except KeyboardInterrupt:
        cleanupPool(pool)


def makeFasta(args, pool=Pool(processes=1)):
    """Finds chimeric sequences from a fasta file and writes them to an accons file.
    :param args: An argparse object with the following parameters:
                    inputFastq	Input Fastq file
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        parallel(runProgramRunner, [ProgramRunner("make.fasta", [args.inputFastq], {"exists": [args.inputFastq]})
                                    ], pool)

    except KeyboardInterrupt:
        cleanupPool(pool)


def prescreen(args, pool=Pool(processes=1)):
    """Prescreens a file for sequences with frameshfits, and logs the frameshifts.

    :param args: An argparse object with the following parameters:
                    aln_out_file        A human-readable output file from vsearch.
                    caln_userout_file   A caln file where each line represents one target/query match.
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        # ensure that our files exist
        helpValidate({"exists": [args.aln_out_file, args.caln_userout_file]})

        # Make the output directory
        makeDir(args.outdir)

        # def screen(vsearch_aln_out_file, vsearch_caln_userout_file):
        screen(args.aln_out_file, args.caln_userout_file)

        # move files
        flagged_seqs_file = "%s.bad" % args.aln_out_file
        move_to_dir(flagged_seqs_file, args.outdir)
    except KeyboardInterrupt:
        cleanupPool(pool)


# Orphaned code
# ========================================================================================================
# ========================================================================================================
# ========================================================================================================


def screenSeqs(args, pool=Pool(processes=1)):
    """Identifies sequences that don't meet specified requirements and writes them to a .accons file.

    :param args: An argparse object with the following parameters:
                    input           Input fasta file to screen.
                    outdir          Directory to dump the output files.
                        ..and any of these filter options:
                    start	        Maximum allowable sequence starting index.
                    end	            Minimum allowable sequence ending index.
                    minlength	    Minimum allowable sequence length.
                    maxlength	    Maximum allowable sequence length.
                    maxambig	    Maxmimum number of allowed ambiguities.
                    maxn	        Maximum number of allowed N's.
                    maxhomop	    Maximum allowable homopolymer length.
                    groups	        Groups file to update.
                    names	        Names file to update.
                    alnReport	    Alignment report to update.
                    contigsReport	Contigs report to update.
                    summaryFile	    SummaryFile to update.

    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # "screen.seqs": "mothur \'#screen.seqs(fasta=%s, %s)\'",

    # identify the input file type
    # TODO Not sure if mothur actually supports these different imput formats for this command.  Documentation is vague.
    try:
        optionString = mothur_buildOptionString(args, mustFilter=True)
        parallel(runProgramRunner, [ProgramRunner("screen.seqs", [args.input, optionString],
                                                  {"exists": [args.input]}, pool)
                                    ])
    except KeyboardInterrupt:
        cleanupPool(pool)


def dropShort(args, poo=Pool(processes=1)):
    good_seqs = []
    for seq in SeqIO.parse(args.inputFasta, "fasta"):
        if len(seq.seq) >= int(args.minLenght):
            good_seqs.append(seq)
        else:
            print "seq %s too short (%s bases)" % (seq.id, len(seq.seq))


# TODO multiproc
def fastxRename(args, pool=Pool(processes=1)):
    """Renames sequences in a fastq file as 1,2,3,...
    :param args: An argparse object with the following parameters:
                    input_f     Forward Fastq Reads
                    input_r     Reverse Fastq Reads
                    outdir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        # Make the output directory, or abort if it already exists
        makeDir(args.outdir)
        printVerbose("\tRenaming sequences")
        # "~/programs/fastx/bin/fastx_renamer -n COUNT -i %s %s"
        rename_outFile_f = os.path.join(args.outdir, os.path.basename(args.input_f) + "_renamed")
        rename_outFile_r = os.path.join(args.outdir, os.path.basename(args.input_r) + "_renamed")
        parallel(runProgramRunner, [
            ProgramRunner("fastx_renamer", [args.input_f, rename_outFile_f], {"exists": [args.input_f]}),
            ProgramRunner("fastx_renamer", [args.input_r, rename_outFile_r], {"exists": [args.input_r]}),
        ], pool)
        printVerbose("\tRenamed %s sequences")
    except KeyboardInterrupt:
        cleanupPool(pool)

def findChimeras(args, pool=Pool(processes=1)):
    """Finds chimeric sequences in a fasta file and writes them to a .accons file.
    :param args: An argparse object with the following parameters:
                    inputFasta	Cleaned inputs File
                    program     Program for detecting and removing chimeras. Default is uchime
                    ...and exactly one of the following:
                    namesFile	.names file to update. See <http://www.mothur.org/wiki/Name_file>
                    refDB       database file to update.
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        supported_programs = ["uchime"]
        if args.program in supported_programs:
            makeDir(args.outdir)
            if args.program == "uchime":

                # find the chimeras (but don't remove them yet)
                # "chmimera.uchime": "mothur #chimera.uchime(fasta=%s[, name=%s] )"
                refs = []
                ref_type = "database"
                if args.refdb is None:
                    ref_type = "name"
                    refs = getInputs(args.names, "*.names")
                else:
                    refs = getInputs(args.refdb)

                inputs = getInputs(args.input, "*_seeds.fasta")

                # expect a 1:1 or 1:n mapping from reference files to input files
                if len(refs) != len(inputs) or len(refs) != 1:
                    print "Error: Did not find correct number of reference files."
                    print "You must provide either exactly one reference file for all inputs,"
                    print "or one reference file per input"
                    exit()

                # if a 1:n mapping from reference files to input files, replicate the reference file name
                if len(refs) == 1:
                    ref_file = refs[0]
                    refs = [ref_file] * len(inputs)

                arg_strings = zip(inputs, refs)

                # mothur "#chimera.uchime(fasta=seeds.fasta, name=seeds.names)"
                parallel(runProgramRunner,
                         [ProgramRunner("chmimera.uchime", [input_, ref_type, ref], {"exists": [input_]})
                          for input_, ref in arg_strings], pool)

            else:
                raise Exception("Unknown program %s for chimera detection or removal" % args.program)
            # Remove from the accon sequences from the input file
            accnos = ["%s/%s.denovo.uchime.accnos" % (getDir(input_), getFileName(input_)) for input_ in inputs]

            # Remove from the accon sequences from the input file
            # removeChimeras.py  seeds.fasta seeds.uchime.accnos
            parallel(runProgramRunner, [ProgramRunner("remove.seqs", [accnos_file, "fasta",
                                                                      "%s.fasta" % ".".join(
                                                                          accnos_file.split(".")[:-3])],
                                                      {"exists": [input_, accnos_file]})
                                        for accnos_file in accnos
                                        if os.path.getsize(accnos_file)], pool)

        out_files = []
        out_file_patts = ["*.accnos", "*.chimeras", "*.pick.fasta"]
        for file_patt in out_file_patts:
            out_files = getInputs(os.path.dirname(args.input), file_patt, critical=False)
        out_files += getInputs(".", "mothur.*.logfile", critical=False)
        bulk_move_to_dir(out_files, args.outdir)

    except KeyboardInterrupt:
        cleanupPool(pool)


def removeSeqs(args, pool=Pool(processes=1)):
    """Removes specific sequences (in an .accons file) from an input file.
    :param args: An argparse object with the following parameters:
                    accnosFile  List of sequence names to remove
                    outdir      Directory to put the output files
                    ...and one of the following input files to clean
                    fasta	    Fasta or fastq file to be cleaned
                    list	    Fasta file to be cleaned.  See <http://www.mothur.org/wiki/List_file>
                    groups  	Fasta file to be cleaned.  See <http://www.mothur.org/wiki/Group_file>
                    names       Fasta file to be cleaned.  See <http://www.mothur.org/wiki/Name_file>
                    count	    Fasta file to be cleaned.  See <http://www.mothur.org/wiki/Count_File>
                    alnReport  	Fasta file to be cleaned.  See <http://www.mothur.org/wiki/Remove.seqs#alignreport_option>
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # "chmimera.uchime":  "mothur \'#remove.seqs(accnos=%s, (fasta=%s|list=%s|groups=%s|names=%s|count=%s|
    #                                                           alignreport=%s)\'",
    # identify the input file type
    supported_file_types = ['fasta', 'list', 'groups', 'names', 'count', 'alnReport']
    if args.filetype not in supported_file_types:
        print "Error: Unsupported input filetype.  Supported types are:"
        for file_type in supported_file_types:
            print file_type
        exit()

    try:
        makeDir(args.outdir)
        accnos = getInputs(args.accnos, "*.accnos", ignoreEmptyFiles=False)
        inputs = getInputs(args.input, "*.%s" % args.filetype, ignoreEmptyFiles=False)
        if len(inputs) != len(accnos):
            print "ERROR: The number of input and .ACCNOS files do not match."
            exit()

        # Remove from the accon sequences from the input file
        zipped_pairs_ = zip(inputs, accnos)
        parallel(runProgramRunner, [ProgramRunner("remove.seqs", [accnos_, args.filetype, input_],
                                                  {"exists": [input_, accnos_]}) for input_, accnos_ in zipped_pairs_],
                 pool)

        # Get the output file name with no directory prefix
        output_files = getInputs("*.pick")
        bulk_move_to_dir(output_files, args.outdir)

    except KeyboardInterrupt:
        cleanupPool(pool)


def align_macse(args, pool=Pool(processes=1)):
    """Aligns sequences by iteratively adding them to a known good alignment.

     :param args: An argparse object with the following parameters:
                    db                  Database against which to align and filter reads
                    samplesDir          Directory containig the samples to be cleaned
                    outdir              Directory where outputs will be saved
     :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        makeDir(args.outdir)
        inputs = getInputs(args.input)

        logging.debug("\t %s Aligning reads using MACSE")
        # Aligns sequences by iteratively adding them to a known good alignment.
        #
        # "macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
        #                                    \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
        #                                    -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
        #                                    -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",

        parallel(runProgramRunner, [ProgramRunner("macse_align",
                                                  [args.db, args.db, input_] +
                                                  ["%s/%s" % (args.outdir, getFileName(input_))] * 3,
                                                  {"exists": [input_, args.db]}) for input_ in inputs], pool)

        logging.debug("\t %s Processing MACSE alignments")
        # Removes non-nucleotide characters in MACSE aligned sequences for all fasta files in the samples directory
        #
        #    (the samplesDir argument).
        # "macse_format":     "java -jar " + programPaths["MACSE"] + "  -prog exportAlignment -align \"%s\" \
        #       -charForRemainingFS - -gc_def 5 -out_AA \"%s\" -out_NT \"%s\" -statFile \"%s\""

        parallel(runProgramRunner, [ProgramRunner("macse_format",
                                                  ["%s/%s_NT" % (args.outdir, getFileName(input_)),
                                                   "%s/%s_AA_macse.fasta" % (args.outdir, getFileName(input_)),
                                                   "%s/%s_NT_macse.fasta" % (args.outdir, getFileName(input_)),
                                                   "%s/%s_macse.csv" % (args.outdir, getFileName(input_))],
                                                  {"exists": [input_]}) for input_ in inputs], pool)

        printVerbose("\tCleaning MACSE alignments")
        # Remove the reference sequences from the MACSE files and remove the non nucleotide characters from the
        #       sequences.
        macse_outputs = ["%s/%s_NT_macse.fasta" % (args.outdir, getFileName(input_)) for input_ in inputs]
        cleaned_output_files = ["%s/MACSE_OUT_.part%d" % (args.outdir, i) for i in range(len(macse_outputs))]
        input_pairs = zip(macse_outputs, cleaned_output_files)
        # TODO removeMacseRefs reads the database with each process, perhpas pass a list of names instead?
        # removeMacseRefs(file_to_clean, reference_file, output_file_name):
        parallel(runPython, [(removeMacseRefs, macse_output, args.db, output_name)
                             for macse_output, output_name in input_pairs], pool)

        # cat the files
        joinFiles(cleaned_output_files, "%s/MACSE_OUT_MERGED.fasta" % args.outdir)

    except KeyboardInterrupt:
        cleanupPool(pool)



# Todo remove this
def test(args, pool=Pool(processes=1)):
    path = ProgramRunner.programPaths[args.program]
    print "PATH: " + path
    subprocess.call(os.path.expanduser(path))
