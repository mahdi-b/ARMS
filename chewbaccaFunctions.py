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
from parsers.parseMHAPOut import findUnmatchedSeqs
from parsers.parseVSearchout import parseVSearchout
from renamers.renameWithReplicantCounts import renameWithReplicantCounts
from renamers.renameSequences import serialRename
from renamers._renameWithCount import renameSequencesWithCount
from renamers.renameWithoutCount import removeCountsFromName, removeCountsFromNamesFile
from renamers.updateNames import updateNames
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
            makeDirOrdie(args.outdir)
            # Run and get a list of output files
            supported_programs[args.program](args, pool)

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
        makeDirOrdie(args.outdir)
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

        printVerbose("\tAssembling reads with pear")
        debugPrintInputInfo(inputs, "assemble")
        parallel(runProgramRunnerInstance, [ProgramRunner("pear", [forwards, reverse,
                                                           "%s/%s_%s" % (args.outdir, args.name, getFileName(forwards)),
                                                                   args.threads],
                                                          {"exists": [forwards, reverse]})
                                            for forwards, reverse in inputs], pool)

        printVerbose("Done assembling sequences...")

        # Grab all the auxillary files (everything not containing ".assembled."
        aux_files = getInputs(args.outdir, "*", "*.assembled.*", ignore_empty_files=False)
        # make aux dir for extraneous files and move them there
        bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    except KeyboardInterrupt:
        cleanupPool(pool)


# TODO handle file cleanup
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
        printVerbose("\tAssembling reads with mothur")
        parallel(runProgramRunnerInstance, [ProgramRunner("trimmomatic", [args.forward, args.reverse, args.bdiffs, args.pdiffs,
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
        makeDirOrdie(args.outdir)

        files_to_split = getInputs(args.input)

        file_id = range(len(files_to_split))
        file_id_pairs = zip(files_to_split, file_id)
        printVerbose("Demuxing sequences...")
        debugPrintInputInfo(files_to_split, "demux")
        parallel(runProgramRunnerInstance, [ProgramRunner("barcode.splitter",
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
        makeDirOrdie(args.outdir)

        # Gather input files
        inputs = getInputs(args.input)
        debugPrintInputInfo(inputs, "rename")

        printVerbose("Renaming sequences...")
        # Run serialRename in parallel
        parallel(runPythonInstance,
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
            makeDirOrdie(args.outdir)
            # Run and get a list of output files
            supported_programs[prog](args, pool)
        else:
            print "Invalid program choice.  Supported programs are " + str(keys)
    except KeyboardInterrupt:
        cleanupPool(pool)


def trim_flexbar(args, pool=Pool(processes=1)):
    # "flexbar":  "flexbar -r \"%s\" -t \"%s\" -ae \"%s\" -a \"%s\"",
    try:
        makeDirOrdie(args.outdir)

        temp_file_name_template = "%s/temp_%s"
        debarcoded_file_name_template = "%s/%s_debarcoded"

        # TODO those exists validators dont really need ot be there since we globbed our files
        # TODO these should write to the input directory, then return the file names, so the parent can move and rename
        printVerbose("Trimming barcodes and adapters with flexbar")
        inputs = getInputs(args.input)
        debugPrintInputInfo(inputs, "trim adapters from")
        # Trim the left
        parallel(runProgramRunnerInstance, [ProgramRunner("flexbar",
                                                          [input_file,
                                                   temp_file_name_template % (args.outdir, strip_ixes(input_file)),
                                                   "LEFT", args.barcodes],
                                                          {"exists": [input_file]}) for input_file in inputs], pool)

        temp_files = getInputs(args.outdir, "temp_*")
        debugPrintInputInfo(temp_files, "debarcode")

        # Trim the right
        parallel(runProgramRunnerInstance, [ProgramRunner("flexbar",
                                                          [input_file,
                                                   debarcoded_file_name_template % (
                                                       args.outdir, strip_ixes(input_file)[5:]),
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
        debugPrintInputInfo(inputs, "trim")
        parallel(runProgramRunnerInstance, [ProgramRunner("trim.seqs", [input_, args.oligos],
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
        makeDirOrdie(args.outdir)

        printVerbose("Cleaning sequences with Trimmomatic...")
        inputs = getInputs(args.input)
        debugPrintInputInfo(inputs, "clean")
        parallel(runProgramRunnerInstance, [ProgramRunner("trimmomatic", [input_,
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
        makeDirOrdie(args.outdir)
        # TODO debug/verbose for conversion
        # grab all the _cleaned files and turn the fastqs into fastas
        input_fqs = getInputs(args.input)
        output_fas = []
        fastq_exts = {".fq", ".fastq"}

        # Generate a list of files to translate from fastq to fasta
        for input_ in input_fqs:
            if os.path.splitext(input_)[1].lower() in fastq_exts:
                base_name = os.path.basename(input_)
                file_name = os.path.splitext(base_name)[0]
                fasta_name = "%s/%s.fa" % (args.outdir, file_name)
                translateFastqToFasta(input_, fasta_name)
                output_fas.append(fasta_name)

        # zip those names to parallel array s
        file_args = zip(input_fqs, output_fas)
        # Translate fastq files in parallel
        parallel(runPythonInstance, [(translateFastqToFasta, input_fq, output_fa) for input_fq, output_fa in file_args], pool)

        debugPrintInputInfo(output_fas, "searched for replicants")
        printVerbose("Searching for identical reads...")
        # parse the number of identical reads included in each sequence and write them to the
        #       {sample_file}_parsed.out
        # ls *cleaned.fasta  | parallel  "~/ARMS/bin/usearch7.0.1090 -derep_fulllength {} -output {/.}_derep.fa
        #                                                   -uc {/.}_uc.out"
        parallel(runProgramRunnerInstance, [ProgramRunner("usearch", [input_,
                                                              "%s/%s_derep.fa" % (
                                                                  args.outdir, strip_ixes(getFileName(input_))),
                                                              "%s/%s_uc.out" % (
                                                                  args.outdir, strip_ixes(getFileName(input_)))],
                                                          {"exists": [args.outdir, input_]})
                                            for input_ in output_fas], pool)
        printVerbose("Done searching.")

        # Collect .uc files
        input_ucs = getInputs(args.outdir, "*_uc.out")
        debugPrintInputInfo(input_ucs, "parsed")
        printVerbose("Parsing search logs for seed sequences...")
        # Collect seeds
        # #ls *uc.out | parallel "python ~/ARMS/bin/getSeedSequences.py {} {.}_parsed.out"
        parallel(runPythonInstance,
                 [(getSeedSequences, input_, "%s/%s.names" % (args.outdir, strip_ixes(getFileName(input_)))) for
                  input_ in input_ucs], pool)
        printVerbose("Done parsing search logs.")

        # Rename the sequences to include the the number of identical reads.
        #  Ex. in the fasta file {sample_file}_derep_renamed.fa, read 123_10 indicate that for sequence
        # which id is 123, there 10 sequences that identical to it and which were discarded.
        # ls * cleaned.fasta | parallel "python ~/ARMS/bin/_renameWithCount.py
        #                   {/.}_derep.fa {/.}_uc_parsed.out {/.}_derep_renamed.fa"
        input_fa = getInputs(args.outdir, "*_derep.fa")
        names_file = getInputs(args.outdir, "*.names")
        inputs = zip(input_fa, names_file)
        debugPrintInputInfo(inputs, "rename, names file")

        # renameSequencesWithCount(input_fasta, count_file, outfile):
        printVerbose("Renaming sequences with replication counts...")
        parallel(runPythonInstance,
                 [(renameSequencesWithCount,
                   fasta_file, count_file,
                   "%s/%s_derepCount.fa" % (args.outdir, strip_ixes(getFileName(fasta_file))))
                  for fasta_file, count_file in inputs], pool)
        printVerbose("Done renaming sequences.")

        aux_dir = makeAuxDir(args.outdir)
        aux_files = getInputs(args.outdir, "*", "*_derepCount.*")
        bulk_move_to_dir(aux_files, aux_dir)
        # separate names folder
        bulk_move_to_dir(getInputs(aux_dir, "*.names"), makeDirOrdie("%s_%s" % (args.outdir, "names_files")))

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

        makeDirOrdie(args.outdir)

        inputs = getInputs(args.input)
        debugPrintInputInfo(inputs, "aligned")

        printVerbose("Aligning sequences with mothur against BIOCODE")
        # ls *debarcoded_cleaned_derep_renamed.fa | parallel  \
        #   mothur '"#align.seqs(candidate={}, template=~/ARMS/data/BIOCODETEMPLATE, flip=t)"'
        #
        # "align.seqs": "mothur \'#align.seqs(candidate=%s, template=%s, flip=t)\'",
        parallel(runProgramRunnerInstance, [ProgramRunner("align.seqs", [input_, args.ref],
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
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # def splitK(inputFasta, prefix, nbSeqsPerFile, filetype):
    try:
        # Make the output directory, or abort if it already exists
        makeDirOrdie(args.outdir)

        # Gather input files
        inputs = getInputs(args.input)

        # Run renamer in parallel
        parallel(runPythonInstance, [
            (splitK, input_, "%s/%s" % (args.outdir, strip_ixes(getFileName(input_))), args.chunksize, args.filetype)
            for input_ in inputs], pool)

    except KeyboardInterrupt:
        cleanupPool(pool)


def merge(args, pool=Pool(processes=1)):
    """Blindly concatenates files in a directory.
       :param args: An argparse object with the following parameters:
                       input        Cleaned inputs File.
                       outdir       Directory where outputs will be saved.
                       name         Name prefix for the merged file.
                       fileext      Output file extension.  e.g 'fasta', 'fastq', 'txt'
       :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
       """
    try:
        makeDirOrdie(args.outdir)
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
        makeDirOrdie(args.outdir)
        inputs = getInputs(args.input)

        # ungap(file_to_clean, output_file_name, gap_char, file_type):
        parallel(runPythonInstance, [(ungap, input_, "%s/%s_cleaned.%s" % (args.outdir, strip_ixes(input_), 'fasta'),
                                      args.gapchar, args.fileext)
                                     for input_ in inputs], pool)
    except KeyboardInterrupt:
        cleanupPool(pool)


def cluster(args, pool=Pool(processes=1)):
    # TODO what is the correct default behavior if no names file is supplied? Currently singletons.
    """Clusters sequences.
    :param args: An argparse object with the following parameters:
                    input       A file or folder containing fasta files to cluster.
                    output      The output directory results will be written to.
                    names       A names file or folder containing names files that describe the input.
                                Note: if no names file is supplied, then entries in the fasta file are assumed to be
                                    singleton sequences.

                    namesFile	Reference .names file. See <http://www.mothur.org/wiki/Name_file>
                    stripcounts If False, leaves any abbundance suffixes (_###) intact.  This may hurt clustering.
                                    Defaults to True.
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        makeDirOrdie(args.outdir)
        inputs = getInputs(args.input)

        # REMOVES COUNTS FROM SEQUENCE NAMES IN ORDER TO CLUSTER PROPERLY
        # strip counts if we need to.
        if args.stripcounts:
            printVerbose("Removing counts from sequence names...")
            debugPrintInputInfo(inputs, "renamed")
            parallel(runPythonInstance,
                     [(removeCountsFromName, input_, "%s/%s_uncount.fasta" % (args.outdir, strip_ixes(input_)), 'fasta')
                      for input_ in inputs], pool)
            printVerbose("Done removing counts.")

            # Grab the cleaned files as input
            inputs = getInputs(args.outdir, "*_uncount.fasta")

        # DEREPLICATE ONE MORE TIME
        printVerbose("Dereplicating before clustering...")
        debugPrintInputInfo(inputs, "dereplicated")
        parallel(runProgramRunnerInstance, [ProgramRunner("vsearch.derep",
                                                          [input_, "%s/%s_derep.fasta" % (args.outdir, strip_ixes(input_)),
                                                   "%s/%s_uc.out" % (args.outdir, strip_ixes(input_))],
                                                          {"exists": [input_]}) for input_ in inputs], pool)
        printVerbose("Done dereplicating")
        # generates a .names file named _uc_parsed.out
        # python getSeedSequences.py uc.out uc_parsed.out
        input_ucs = getInputs(args.outdir, "*_uc.out")
        printVerbose("Generating a names file from dereplication.")
        debugPrintInputInfo(inputs, "parsed (into a names file)")
        parallel(runPythonInstance,
                 [(getSeedSequences, input_, "%s/%s_derep.names" % (args.outdir, strip_ixes(input_)))
                  for input_ in input_ucs], pool)

        most_recent_names_files = getInputs(args.outdir, "*_derep.names")

        # UPDATE THE NAMES FILES
        if args.namesfile is not None:
            # Grab the old names file and the dereplicated names file
            old_names_files = getInputs(args.namesfile)
            derep_names_files = getInputs(args.outdir, "*_derep.names")

            printVerbose("Updating .names files with dereplicated data")
            logging.debug("%d Reference (old) names files to be read:" % len(old_names_files))
            logging.debug(str(old_names_files))
            logging.debug("%d Dereplicated (new) names files to be read:" % len(derep_names_files))
            logging.debug(str(derep_names_files))
            # updateNames (old_names_files, new_names_files, updated)
            updateNames(old_names_files, derep_names_files, args.outdir, "precluster")
            most_recent_names_files = getInputs(args.outdir, "precluster*")
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

        # CLUSTER
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
        inputs = getInputs(args.outdir, "*_counts.fasta")
        debugPrintInputInfo(inputs, "clustered")
        parallel(runProgramRunnerInstance, [ProgramRunner("swarm", [input_,
                                                            "%s/%s_clustered.names" % (args.outdir, strip_ixes(input_)),
                                                            "%s/%s_clustered_uclust" % (
                                                                args.outdir, strip_ixes(input_)),
                                                            "%s/%s_clustered_seeds" % (
                                                                args.outdir, strip_ixes(input_))],
                                                          {"exists": [input_]}) for input_ in inputs], pool)

        # REMOVE COUNTS FROM CLUSTERING NAMES FILE
        # Grab the current names file and the new clustered names file (which needs to be cleaned)
        clustered_names_files = getInputs(args.outdir, "*_clustered.names")
        # Remove counts from the clustering names files
        printVerbose("Cleaning the .names file from clustering")
        debugPrintInputInfo(clustered_names_files, "cleaned")
        parallel(runPythonInstance,
                 [(removeCountsFromNamesFile, input_, "%s/%s_uncount.names" % (args.outdir, strip_ixes(input_)))
                  for input_ in clustered_names_files], pool)

        cleaned_clustered_names_files = getInputs(args.outdir, "*clustered_uncount.names")
        # UPDATE THE NAMES FILES WITH NEW CLUSTERS
        printVerbose("Updating .names files with clustering data")
        # updateNames (old_names_files, new_names_files, updated)
        print most_recent_names_files
        print cleaned_clustered_names_files
        updateNames(most_recent_names_files, cleaned_clustered_names_files, args.outdir, "postcluster")
        printVerbose("Done updating .names files.")

        # Convert the seeds files to uppercase (swarm writes in lowercase)
        inputs = getInputs(args.outdir, "*_seeds")
        parallel(runPythonInstance,
                 [(capitalizeSeqs, input_, "%s.fasta" % input_) for input_ in inputs], pool)
        # delete seeds file
        for input_ in inputs:
            os.remove(input_)
            inputs = getInputs(args.outdir, "*_seeds.fasta")
            parallel(runPythonInstance,
                     [(seedToNames, input_, "%s/%s.names" % (args.outdir, getFileName(input_))) for input_ in inputs],
                     pool)

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
        db_string = os.path.expanduser("~/ARMS/data/BiocodePASSED_SAP.txt")
        aln_user_string = ""
        makeDirOrdie(args.outdir)

        # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
        #       --userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt

        # expecting a fasta to annotate
        inputs = getInputs(args.input)
        # input, db, userout, alnout, aln_userfield
        parallel(runProgramRunnerInstance, [ProgramRunner("vsearch.usearch_global",
                                                          [input_, db_string, "%s/%s.out" % (args.outdir, strip_ixes(input_)),
                                                   "%s/%s.alnout" % (args.outdir, strip_ixes(input_)), aln_user_string],
                                                          {"exists": [input_]}) for input_ in inputs], pool)
        # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
        # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
        #
        # parseVSearchout.py ~/ARMS/data/BiocodePASSED_SAP_tax_info.txt vsearch_out  parsed_BIOCODE.out 97 85
        parallel(runPythonInstance, [(parseVSearchout, os.path.expanduser("~/ARMS/data/BiocodePASSED_SAP_tax_info.txt"),
                              "%s/%s.out" % (args.outdir, strip_ixes(input_)),
                              "%s/%s_result.out" % (args.outdir, strip_ixes(input_)),
                                      97, 85) for input_ in inputs], pool)
        # Gather and move auxillary files
        aux_files = getInputs(args.outdir, "*", "*_result.out")
        bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

    except KeyboardInterrupt:
        cleanupPool(pool)


def queryNCBI(args, pool=Pool(processes=1)):
    """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

    :param args: An argparse object with the following parameters:
                    input       Input file/folder with fasta sequences
                    outdir      Directory to put the output files
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        coi_db_string = os.path.expanduser("~/ARMS/refs/COI.fasta")
        ncbi_db_string = os.path.expanduser("~/ARMS/refs/ncbi.db")
        aln_user_string = "--userfields query+target+id+alnlen+qcov"
        makeDirOrdie(args.outdir)

        # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
        #		--userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt

        # expecting a fasta to annotate
        inputs = getInputs(args.input)
        # input, db, userout, alnout, aln_userfield
        parallel(runProgramRunnerInstance, [ProgramRunner("vsearch.usearch_global",
                                                          [input_, coi_db_string,
                                                   "%s/%s.out" % (args.outdir, strip_ixes(input_)),
                                                   "%s/%s.alnout" % (args.outdir, strip_ixes(input_)), aln_user_string],
                                                          {"exists": [input_]}) for input_ in inputs], pool)
        # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
        # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
        #
        # parseVSearchToTaxa(vsearch_out, ncbi_db, min_coverage, min_similarity)> parsed_nt.out
        parallel(runPythonInstance, [(parseVSearchToTaxa, "%s/%s.out" % (args.outdir, strip_ixes(input_)), ncbi_db_string,
                              "%s/%s_result.out" % (args.outdir, strip_ixes(input_)),
                                      97, 85) for input_ in inputs], pool)

        # Gather and move auxillary files
        aux_files = getInputs(args.outdir, "*", "*_result.out")
        bulk_move_to_dir(aux_files, makeAuxDir(args.outdir))

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
        parallel(runProgramRunnerInstance, [ProgramRunner("make.fastq", [args.inputFasta, args.inputQual],
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
        parallel(runProgramRunnerInstance, [ProgramRunner("make.fasta", [args.inputFastq], {"exists": [args.inputFastq]})
                                            ], pool)

    except KeyboardInterrupt:
        cleanupPool(pool)


def prescreen(args, pool=Pool(processes=1)):
    """Prescreens a file for sequences with frameshfits, and logs the frameshifts.

    :param args: An argparse object with the following parameters:
    input           input file/folder with fasta sequences
                    outdir      Directory to put the output files
                    aln_out_file        A human-readable output file from vsearch.
                    caln_userout_file   A caln file where each line represents one target/query match.
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        # ensure that our files exist
        helpValidate({"exists": [args.aln_out_file, args.caln_userout_file]})

        # Make the output directory
        makeDirOrdie(args.outdir)

        # def screen(vsearch_aln_out_file, vsearch_caln_userout_file):
        screen(args.aln_out_file, args.caln_userout_file)

        # move files
        flagged_seqs_file = "%s.bad" % args.aln_out_file
        move_to_dir(flagged_seqs_file, args.outdir)
    except KeyboardInterrupt:
        cleanupPool(pool)

#TODO doc and debug
def minhash(args, pool=Pool(processes=1)):
    """Queries a compiled database using minhashes to determine the orientation of a file.

    :param args: An argparse object with the following parameters:
                    input       File/folder with fasta sequences
                    outdir      Directory to put the output files
                    rebuilddb   If true, recompile the DB
                    dbfasta     The fasta file listing reference sequences (the DB to compile)
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    #match once, then match again with stragglers

    database_dir = os.path.expanduser("~/ARMS/data/")
    makeDirOrdie(database_dir, orDie=False)
    compiled_data_base = getInputs("%s%s.dat" % (database_dir, getFileName(args.dbfasta)), critical=False)

    # if the db doesnt exist, or forced rebuild,
    if len(compiled_data_base) == 0 or args.rebuilddb:
        print "Compiling database..."
        # (Re)Build the Database from a fasta
        db_source_fastas = getInputs(args.dbfasta)
        debugPrintInputInfo(db_source_fastas, "be compiled into a DB.")
        if len(db_source_fastas) == 0:
            print("Error: Database source fasta(s) not found!")
            sys.exit()

        #"min.hash.build.db": "java -jar ~/ARMS/programs/mhap/mhap-2.1.jar --store-full-id -p \"%s\" -q \"%s\"",
        parallel(runProgramRunnerInstance, [ProgramRunner("min.hash.build.db",
                                                          [source, database_dir],
                                                          {"exists": []})
                                            for source in db_source_fastas], pool)
    makeDirOrdie(args.outdir)
    # Convert the sensativity list to a reverse-sorted list of integers
    sensativity_list = map(int, args.sensativities.split(","))
    sensativity_list.sort(reverse=True)

    query_fastas = getInputs(args.input)
    debugPrintInputInfo(query_fastas, "be queried against the db.")
    for sensitivity in sensativity_list:

        #"min.hash.search": "java -jar ~/ARMS/programs/mhap/mhap-2.1.jar --store-full-id -s \"%s\" -q \"%s\" \
        #                                       --no-self --num-min-matches %d > \"%s\"",
        printVerbose("Querying DB with sensativity %d" % sensitivity)
        parallel(runProgramRunnerInstance, [ProgramRunner("min.hash.query",
                                              [compiled_data_base[0], query, sensitivity, "%s/%s_%d.out" %
                                                (args.outdir, clip_count(strip_ixes(query)), sensitivity)], {"exists": []})
                                            for query in query_fastas], pool)
        printVerbose("Done with queries.")

        # Diff the found items from the query file,
        mhap_outputs = getInputs(args.outdir, "*_%d.out" % sensitivity, ignore_empty_files=False)

        # def findUnmatchedSeqs(fasta_to_clean, mhap_outfile, unmatched_fasta):
        fasta_mhap_pairs = zip(query_fastas, mhap_outputs)
        debugPrintInputInfo(fasta_mhap_pairs, "parse and screen from the query fastas.")

        printVerbose("Removing identified sequences from query fasta")
        parallel(runPythonInstance, [(findUnmatchedSeqs, fasta_query, mhap_output,
                                      "%s/%s_%d.fasta" % (args.outdir, clip_count(strip_ixes(fasta_query)),sensitivity))
                                     for fasta_query, mhap_output in fasta_mhap_pairs], pool)
        printVerbose("Done removing sequences")

        # update counters for next iteration
        query_fastas = getInputs(args.outdir, "*_%d.fasta" % sensitivity, ignore_empty_files=False)


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
        option_string = mothur_buildOptionString(args, mustFilter=True)
        parallel(runProgramRunnerInstance, [ProgramRunner("screen.seqs", [args.input, option_string],
                                                          {"exists": [args.input]}, pool)
                                            ])
    except KeyboardInterrupt:
        cleanupPool(pool)

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
        makeDirOrdie(args.outdir)
        printVerbose("\tRenaming sequences")
        # "~/programs/fastx/bin/fastx_renamer -n COUNT -i %s %s"
        rename_out_file_f = os.path.join(args.outdir, os.path.basename(args.input_f) + "_renamed")
        rename_out_file_r = os.path.join(args.outdir, os.path.basename(args.input_r) + "_renamed")
        parallel(runProgramRunnerInstance, [
            ProgramRunner("fastx_renamer", [args.input_f, rename_out_file_f], {"exists": [args.input_f]}),
            ProgramRunner("fastx_renamer", [args.input_r, rename_out_file_r], {"exists": [args.input_r]}),
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
            makeDirOrdie(args.outdir)
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
                parallel(runProgramRunnerInstance,
                         [ProgramRunner("chmimera.uchime", [input_, ref_type, ref], {"exists": [input_]})
                          for input_, ref in arg_strings], pool)

            else:
                raise Exception("Unknown program %s for chimera detection or removal" % args.program)
            # Remove from the accon sequences from the input file
            accnos = ["%s/%s.denovo.uchime.accnos" % (getDir(input_), getFileName(input_)) for input_ in inputs]

            # Remove from the accon sequences from the input file
            # removeChimeras.py  seeds.fasta seeds.uchime.accnos
            parallel(runProgramRunnerInstance, [ProgramRunner("remove.seqs", [accnos_file, "fasta",
                                                                      "%s.fasta" % ".".join(
                                                                          accnos_file.split(".")[:-3])],
                                                              {"exists": [accnos_file]})
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
                    list	    List file to be cleaned.  See <http://www.mothur.org/wiki/List_file>
                    groups  	Group file to be cleaned.  See <http://www.mothur.org/wiki/Group_file>
                    names       Names file to be cleaned.  See <http://www.mothur.org/wiki/Name_file>
                    count	    Count file to be cleaned.  See <http://www.mothur.org/wiki/Count_File>
                    alnReport  	Aln file to be cleaned.  See <http://www.mothur.org/wiki/Remove.seqs#alignreport_option>
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
        makeDirOrdie(args.outdir)
        accnos = getInputs(args.accnos, "*.accnos", ignore_empty_files=False)
        inputs = getInputs(args.input, "*.%s" % args.filetype, ignore_empty_files=False)
        if len(inputs) != len(accnos):
            print "ERROR: The number of input and .ACCNOS files do not match."
            exit()

        # Remove from the accon sequences from the input file
        zipped_pairs_ = zip(inputs, accnos)
        parallel(runProgramRunnerInstance, [ProgramRunner("remove.seqs", [accnos_, args.filetype, input_],
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
                    input               Directory containig the samples to be cleaned
                    outdir              Directory where outputs will be saved
     :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        makeDirOrdie(args.outdir)
        inputs = getInputs(args.input)

        logging.debug("\t %s Aligning reads using MACSE")
        # Aligns sequences by iteratively adding them to a known good alignment.
        #
        # "macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
        #                                    \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
        #                                    -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
        #                                    -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",

        parallel(runProgramRunnerInstance, [ProgramRunner("macse_align",
                                                          [args.db, args.db, input_] +
                                                          ["%s/%s" % (args.outdir, getFileName(input_))] * 3,
                                                          {"exists": [input_, args.db]}) for input_ in inputs], pool)

        logging.debug("\t %s Processing MACSE alignments")
        # Removes non-nucleotide characters in MACSE aligned sequences for all fasta files in the samples directory
        #
        #    (the samplesDir argument).
        # "macse_format":     "java -jar " + programPaths["MACSE"] + "  -prog exportAlignment -align \"%s\" \
        #       -charForRemainingFS - -gc_def 5 -out_AA \"%s\" -out_NT \"%s\" -statFile \"%s\""

        parallel(runProgramRunnerInstance, [ProgramRunner("macse_format",
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
        parallel(runPythonInstance, [(removeMacseRefs, macse_output, args.db, output_name)
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
