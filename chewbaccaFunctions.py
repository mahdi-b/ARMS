import  subprocess
import time
from Bio.Seq import Seq
from multiprocessing import Pool
from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner
from converters.fastqToFasta import translateFastqToFasta
from getSeedSequences import getSeedSequences
from renamers.renameSequences import serialRename
from renamers.renameWithCount import renameSequencesWithCount
from utils.splitKperFasta import splitK
from prescreen import screen


def assemble(args, pool=Pool(processes=1)):
    """Assembles reads from two (left and right) fastq files.  Handler for the mothur and pear assemblers.  Validates
        argument dependencies.
    :param args: An argparse object with the following parameters:
                program', type=int, help="The number of threads to use
                name		    Assembled File Prefix
                input_f		    Forward Fastq Reads
                input_r		    Reverse Fastq Reads
                outdir		    Directory where outputs will be saved
                threads         The number of threads to use
                # options for mothur
                oligos          Oligos file with barcode and primer sequences
                bdiffs          # of allowed barcode mismatches
                pdiffs          # of allowed primer mismatches
                # options for pear
                maxlength      Maximum length for assembled sequences
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
            # Relocate all output files
            bulk_move_to_dir(outputs, args.outdir)
        else:
            "Invalid program choice.  Supported programs are " + str(keys)
    except KeyboardInterrupt:
        pool.terminate()


def assemble_pear(args, pool=Pool(processes=1)):
    """Assembles reads from two (left and right) fastq files.
    :param args: An argparse object with the following parameters:
                    name            Textual ID for the data set
                    input_f         Forward Fastq Reads
                    input_r         Reverse Fastq Reads
                    threads         The number of threads to use durring assembly.
                    outdir          Directory where outputs will be saved
                    maxlength      Maximum length for assembled sequences
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """

    # *****************************************************************************************
    # Making the contigs using Pear
    # "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s -m %d"
    try:
        name_error_text = "Forwards reads should include the filename suffix \"_forward.fq\" "\
                    + "or \"R1.fq\".  Reverse reads should include the filename suffix \"_reverse\" or \"R2.fq\"."
        printVerbose("\tAssembling reads with pear")
        threads = 1
        ids = []
        forwards_reads = getInputs(args.input_f, "*_forward.f*")
        forwards_reads += getInputs(args.input_f, "*_R1.f*")
        reverse_reads = getInputs(args.input_r, "*_reverse.f*")
        reverse_reads += getInputs(args.input_r, "*_R2.f*")

        if len(forwards_reads) != len(reverse_reads):
            print "Error: Unequal number of forwards/reverse reads."
            exit()

        if len(forwards_reads) == 0:
            print "Error: No read files found."
            print name_error_text
            exit()

        inputs = zip(forwards_reads, reverse_reads)

        if args.threads is not None and args.threads > 1:
            threads = args.threads

        parallel(runProgramRunner, [ProgramRunner("pear", [forwards, reverse,
                                                           "%s/%s_%s" % (args.outdir, args.name, getFileName(forwards)), threads],
                                                  {"exists": [forwards, reverse]})
                                    for forwards, reverse in inputs], pool)
    except KeyboardInterrupt:
        pool.terminate()
        pool.join()

def assemble_mothur(args, pool=Pool(processes=1)):
    """Finds chimeric sequences from a fasta file and writes them to an accons file.
    :param args: An argparse object with the following parameters:
                    forward	    Forward read fastq file
                    reverse	    Reverse read fastq file
                    oligos	    Oligos file with barcode and primer sequences
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
        return ["%s.trim"]
    except KeyboardInterrupt:
        pool.terminate()


# TODO Document that if samples from site A occur in more than one file, then they will overwrite eachother when you run
# TODO    the barcode splitter.  Seqences from each sample need to be put in their own file.
# TODO fix: if a sample file exists, do sample_1, sample_2 etc..
def splitOnBarcodes(args, pool=Pool(processes=1)):
    """Splits a fasta/fastq file on a set of barcodes.  An output file will be created for each sample, listing all members
        from that sample.
    :param args: An argparse object with the following parameters:
                    inputFile   File to split or all
                    barcodes    Tab delimited files of barcodes and their samples
                    outdir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # TODO MAHDI  Comment the rest of the program to replace the pipeline with
    # 1- split_libraries_fastq.py
    # 2- usearch -fastq_filter
    # "barcode.splitter": "cat \"%s\" | " + programPaths["FASTX"] + 'fastx_barcode_splitter.pl  --bcfile "%s" \
    #                                    -prefix "%s" --suffix %s --bol --mismatches 1',
    try:
        makeDir(args.outdir)
        files_to_split = getInputs(args.input, "*.assembled.f*")
        print files_to_split
        file_id = range(len(files_to_split))
        file_id_pairs = zip(files_to_split, file_id)
        printVerbose("Splitting based on barcodes")
        parallel(runProgramRunner, [ProgramRunner("barcode.splitter",
                                                  [input_, args.barcodes, "%s/" % args.outdir,
                                                   "_splitOut_%d.fastq" % id_],
                                                  {"exists": [input_]}) for input_, id_ in file_id_pairs], pool)

        # wait for output files to be written
        pool.close()
        pool.join()
        printVerbose("Demuxed sequences.")

        # gatehr output files and move them to their final destination
        output_files = enumerateDir(".", "*_splitOut_")
        bulk_move_to_dir(output_files, args.outdir)

    except KeyboardInterrupt:
        pool.terminate()


def trim(args, pool=Pool(processes=1)):
    """Trims the adapter and barcode from each sequence.
    """
    try:
        supported_programs = {
            "flexbar":trim_flexbar,
            "mothur":trim_mothur}
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
            outputs = supported_programs[prog](args, pool)
            # Relocate all output files
            # bulk_move_to_dir(outputs, args.outdir)
        else:
            print "Invalid program choice.  Supported programs are " + str(keys)
    except KeyboardInterrupt:
        pool.terminate()


def trim_flexbar(args, pool=Pool(processes=1)):
    # "flexbar":  "flexbar -r \"%s\" -t \"%s\" -ae \"%s\" -a \"%s\"",
    try:
        makeDir(args.outdir)
        input_file_name_format = "*_renamed.f*"
        input_files = getInputs(args.input, input_file_name_format)
        print("inputfiles:")
        print input_files
        arg = [(input_file, "%s/temp_%s" % (args.outdir, getFileName(input_file)), "LEFT", args.barcodes) for input_file
               in input_files]
        print arg

        printVerbose("Trimming barcodes and adapters with flexbar")
        # TODO those exists validators dont really need ot be there since we globbed our files
        # TODO these should write to the input directory, then return the file names, so the parent can move and rename
        # Trim the left
        parallel(runProgramRunner, [ProgramRunner("flexbar",
                                                  [input_file, "%s/temp_%s" % (args.outdir, getFileName(input_file)),
                                                   "LEFT", args.barcodes],
                                                  {"exists": [input_file]}) for input_file in input_files], pool)

        temp_files = getInputs(args.outdir, "temp_*")
        print input_files
        # Trim the right
        parallel(runProgramRunner, [ProgramRunner("flexbar",
                                                  [input_file,
                                                   "%s/%s_debarcoded" % (args.outdir, getFileName(input_file)[5:]),
                                                   "RIGHT", args.adapters],
                                                  {"exists": [input_file]}) for input_file in temp_files], pool)

        # wait for output files to be written
        pool.close()
        pool.join()
        printVerbose("Demuxed sequences.")

        # gather output files and move them to their final destination
        output_files = enumerateDir(".", "*_debarcoded*")

        return output_files
    except KeyboardInterrupt:
        pool.terminate()


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
        inputs = getInputs(args.input, "*_splitOut_", "*unmatched.*")
        parallel(runProgramRunner, [ProgramRunner("trim.seqs", [input_, args.oligos],
                                                  {"exists": [args.oligos]})
                                    for input_ in inputs], pool)
        printVerbose("Trimmed sequences.")
    except KeyboardInterrupt:
        pool.terminate()


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
        inputs = getInputs(args.input, "*_debarcoded.f*")
        parallel(runProgramRunner, [ProgramRunner("trimmomatic", [input_,
                                                                  "%s/%s_cleaned.fastq" % (
                                                                      args.outdir, getFileName(input_)),
                                                                  args.windowSize, args.quality, args.minLen],
                                                  {"exists": [args.outdir, input_],
                                                   "positive": [args.windowSize, args.quality, args.minLen]})
                                    for input_ in inputs], pool)
    except KeyboardInterrupt:
        pool.terminate()


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
        # grab all the _cleaned files and turn the fastqs into fastas
        input_fqs = getInputs(args.input, "*_cleaned*")
        output_fas = []
        fastq_exts = set([".fq", ".fastq"])

        # Generate a list of files to translate from fastq to fasta
        for input_ in input_fqs:
            if os.path.splitext(input_)[1].lower() in fastq_exts:
                input_dir = os.path.dirname(input_)
                base_name = os.path.basename(input_)
                file_name = os.path.splitext(base_name)[0]
                fasta_name = "%s/%s.fa" % (input_dir, file_name)
                translateFastqToFasta(input_, fasta_name)
                output_fas.append(fasta_name)

        # zip those names to parallel array s
        file_args = zip(input_fqs, output_fas)

        # Translate fastq files in parallel
        parallel(runPython, [(translateFastqToFasta, input_fq, output_fa) for input_fq, output_fa in file_args], pool)

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

        # Collect .uc files
        input_ucs = getInputs(args.outdir, "*_uc.out")

        # Collect seeds
        # #ls *uc.out | parallel "python ~/ARMS/bin/getSeedSequences.py {} {.}_parsed.out"
        parallel(runPython,
                 [(getSeedSequences, input_, "%s/%s_parsed.out" % (args.outdir, strip_ixes(getFileName(input_)))) for
                  input_ in input_ucs], pool)

        # Rename the sequences to include the the number of identical reads.
        #  Ex. in the fasta file {sample_file}_derep_renamed.fa, read 123_10 indicate that for sequence
        # which id is 123, there 10 sequences that identical to it and which were discarded.
        # ls * cleaned.fasta | parallel "python ~/ARMS/bin/renameWithCount.py {/.}_derep.fa {/.}_uc_parsed.out {/.}_derep_renamed.fa"
        input_fa = getInputs(args.input, "*_cleaned.fa")
        input_uc_parsed = getInputs(args.outdir, "*_parsed.out")
        inputs = zip(input_fa, input_uc_parsed)
        # renameSequencesWithCount(input_fasta, count_file, outfile):
        parallel(runPython,
                 [(renameSequencesWithCount,
                   fasta_file, count_file,
                   "%s/%s_derep_renamed.fa" % (args.outdir, strip_ixes(getFileName(fasta_file))))
                  for fasta_file, count_file in inputs], pool)


    except KeyboardInterrupt:
        pool.terminate()


def renameSequences(args, pool=Pool(processes=1)):
    """Renames sequences in a fastq file as 1,2,3,...
    :param args: An argparse object with the following parameters:
                    input     Forward Fastq Reads
                    outdir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        # Make the output directory, or abort if it already exists
        makeDir(args.outdir)

        # Gather input files
        inputs = getInputs(args.input, "*splitOut_*", "unmatched_*")

        # Run renamer in parallel
        parallel(runPython,
                 [(serialRename,
                   input, "%s/%s_renamed%s" % (args.outdir, strip_ixes(getFileName(input)), os.path.splitext(input)[1]),
                   args.filetype)
                  for input in inputs],
                 pool)
        pool.close()
        pool.join()

    except KeyboardInterrupt:
        pool.terminate()


# TODO decide whether to use this function or splitOnBarcodes.  This one uses SeqIO, the above users fastx
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
        inputs = getInputs(args.input, "*.align")

        # Run renamer in parallel
        parallel(runPython, [
            (splitK, input, "%s/%s" % (args.outdir, strip_ixes(getFileName(input))), args.chunksize, args.filetype)
            for input in inputs], pool)
    except KeyboardInterrupt:
        pool.terminate()


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

        inputs = getInputs(args.input, "*_derep.fa")

        # ls *debarcoded_cleaned_derep_renamed.fa | parallel  \
        #   mothur '"#align.seqs(candidate={}, template=~/ARMS/data/BIOCODETEMPLATE, flip=t)"'
        #
        # "align.seqs": "mothur \'#align.seqs(candidate=%s, template=%s, flip=t)\'",
        parallel(runProgramRunner, [ProgramRunner("align.seqs", [input, args.ref],
                                                  {"exists": [args.outdir, args.ref, input]})
                                    for input in inputs], pool)

        out_exts = ["*.align", "*.align.report", "*.accnos"]

        # Wait for files to get written
        pool.close()
        pool.join()

        if len(inputs) > 0:
            indir = os.path.dirname(inputs[0])
            out_files = []
            for ext in out_exts:
                out_files += getInputs(indir, ext)
            bulk_move_to_dir(out_files, args.outdir)

    except KeyboardInterrupt:
        pool.terminate()


def align_macse(args, pool=Pool(processes=1)):
    """Aligns sequences by iteratively adding them to a known good alignment.

     :param args: An argparse object with the following parameters:
                    db                  Database against which to align and filter reads
                    samplesDir          Directory containig the samples to be cleaned
                    outdir              Directory where outputs will be saved
     :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        pool_size = pool._processes
        makeDir(args.outdir)
        inputs = getInputs(args.input)
        printVerbose("\t %s Aligning reads using MACSE")
        #Aligns sequences by iteratively adding them to a known good alignment.
        #
        # "macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
        #                                    \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
        #                                    -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
        #                                    -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",
        done = False
        done = parallel(runProgramRunner, [ProgramRunner("macse_align",
                                                  [args.db, args.db, input_] +
                                                  ["%s/%s" % (args.outdir, getFileName(input_))] * 3
                                                  , {"exists": [input_, args.db]}) for input_ in inputs], pool)
        print("Done writing")
        while not done:
            time.sleep(1000)
        done = False
        printVerbose("\t %s Processing MACSE alignments")
        # Removes non-nucleotide characters in MACSE aligned sequences for all fasta files in the samples directory
        #
        #    (the samplesDir argument).
        # "macse_format":     "java -jar " + programPaths["MACSE"] + "  -prog exportAlignment -align \"%s\" \
        #       -charForRemainingFS - -gc_def 5 -out_AA \"%s\" -out_NT \"%s\" -statFile \"%s\""
        parallel(runProgramRunner, [ProgramRunner("macse_format",
                                                  ["%s/%s_NT" % (args.outdir, getFileName(input)),
                                                   "%s/%s_AA_macse.fasta" % (args.outdir, getFileName(input)),
                                                   "%s/%s_NT_macse.fasta" % (args.outdir, getFileName(input)),
                                                   "%s/%s_macse.csv" % (args.outdir, getFileName(input))],
                                                  {"exists": []}) for input_ in os.listdir(args.inputs)], pool)
        while not done:
            time.sleep(1000)

        pool = Pool(pool_size)
        printVerbose("\tCleaning MACSE alignments")
        # Remove the reference sequences from the MACSE files and remove the non nucleotide characters from the sequences.
        #       we need the datbase seq. names to remove them from the results files
        # TODO:IMPORTANT: Merge the files before doing this.
        good_seqs = []
        dbSeqNames = SeqIO.to_dict(SeqIO.parse(args.db, "fasta")).keys()
        print "Will be processing %s samples " % len(inputs)
        i = 0

        # TODO: write parallel script for this
        for sample in inputs:
            nt_macse_out = "%s/%s_NT_macse.fasta" % (args.outdir, sample)
            for mySeq in SeqIO.parse(nt_macse_out, 'fasta'):
                if mySeq.id not in dbSeqNames:
                    mySeq.seq = Seq(str(mySeq.seq[2:]).replace("-", ""))  # remove the !! from the beginning
                    good_seqs.append(mySeq)
            print "completed %s samples" % i
            i += 1
        SeqIO.write(good_seqs, open(os.path.join(args.outdir, "MACSE_OUT_MERGED.fasta"), 'w'), 'fasta')

        #printVerbose("\t%s sequences cleaned, %s sequences retained, %s sequences discarded" % (1, 1, 1))

        pool.join()
    except KeyboardInterrupt:
        pool.terminate()


def findChimeras(args, pool=Pool(processes=1)):
    """Finds chimeric sequences in a fasta file and writes them to a .accons file.
    :param args: An argparse object with the following parameters:
                    inputFasta	Cleaned inputs File
                    program     Program for detecting and removing chimeras. Default is uchime
                    ...and exactly one of the following:
                    namesFile	Reference .names file. See <http://www.mothur.org/wiki/Name_file>
                    refDB       Reference database file.
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        if args.program == "uchime":
            # find the chimeras (but don't remove them yet)
            # "chmimera.uchime": "mothur #chimera.uchime(fasta=%s, (name=%s | reference=%s) )",
            referenceString = ""
            refFile = ""
            if args.namesFile:
                referenceString = "names=%s" % args.namesFile
                refFile = args.namesFile
            else:
                referenceString = "reference=%s" % args.refDB
                refFile = args.refDB

            parallel(runProgramRunner, [ProgramRunner("chmimera.uchime",
                                                      [args.inputFile, referenceString],
                                                      {"exists": [args.inputFasta, refFile]})], pool)

        else:
            raise Exception("unknown program %s for chimera detection or removal" % args.program)

    except KeyboardInterrupt:
        pool.terminate()


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
    inputFileType = ""
    inputFile = ""
    try:
        # Remove from the accon sequences from the input file
        parallel(runProgramRunner, [ProgramRunner("remove.seqs", [args.accnosFile, args.inputFile],
                                                  {"exists": [args.accnosFile, inputFile]}, pool)
                                    ])
        # Get the output file name with no directory prefix
        inputFileName = os.path.basename(inputFile)
        splitInputFileName = inputFileName.split(".")
        splitInputFileName.insert(-1, "pick")
        pickOutFile = (".").join(splitInputFileName)

        # TODO move output '*.pick.filetype' file to 'outputDir/*.pick.filetype'
        move("%s/%s" % (os.path.dirname(inputFile), pickOutFile),
             "%s/%s" % (args.outdir, pickOutFile))
        printVerbose("\t Removed %s target sequences")

    except KeyboardInterrupt:
        pool.terminate()


def screenSeqs(args, pool=Pool(processes=1)):
    """Identifies sequences that don't meet specified requirements and writes them to a .accons file.

    :param args: An argparse object with the following parameters:
                    inputfile       Input fasta file to screen.
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
        parallel(runProgramRunner, [ProgramRunner("screen.seqs", [args.inputfile, optionString],
                                                  {"exists": [args.inputfile]}, pool)
                                    ])

    except KeyboardInterrupt:
        pool.terminate()


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
        pool.terminate()


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
        pool.terminate()


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
        pool.terminate()


# Orphaned code
# ========================================================================================================
# ========================================================================================================
# ========================================================================================================
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
        pool.terminate()


# Todo remove this
def test(args, pool=Pool(processes=1)):
    path = ProgramRunner.programPaths[args.program]
    print "PATH: " + path
    subprocess.call(os.path.expanduser(path))
