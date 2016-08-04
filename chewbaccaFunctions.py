import glob
import re
from Bio.Seq import Seq
from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner
from multiprocessing import Pool
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
                # options
                # for mothur
                oligos          Oligos file with barcode and primer sequences
                bdiffs          # of allowed barcode mismatches
                pdiffs          # of allowed primer mismatches
                # for pear
                maxlength      Maximum length for assembled sequences
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    supportedPrograms = {
        "mothur": assemble_mothur,
        "pear": assemble_pear
    }
    keys = supportedPrograms.keys()
    if args.program in keys:
        # Validate argument dependencies
        if args.program == "mothur":
            if not args.oligos:
                "-g parameter required for mothur assembler"
        # Make the output directory
        makeDir(args.outdir)
        # Run and get a list of output files
        outputs = supportedPrograms[args.program](args, pool)
        # Relocate all output files
        bulk_move_to_dir(outputs, args.outdir)
    else:
        "Invalid program choice.  Supported programs are " + keys

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
        printVerbose("\tAssembling reads with pear")
        # files must be named name_forward or name_reverse
        forwardPrefix = "_forward"
        reversePrevix = "_reverse"
        names = []

        forwards_reads = getInputs(args.input_f)
        reverse_reads = getInputs(args.input_r)

        if len(forwards_reads) != len(reverse_reads):
            print "Error: Unequal number of forwards/reverse reads."
            exit()

        if len(forwards_reads) == 0:
            print "Error: No read files found"
            exit()

        regex_matches = [re.search(r'(.*)_forward\..*', filename.split("/")[-1]) for filename in forwards_reads]
        for match in regex_matches:
            if match != None:
                print(match)
                names.append(match.groups(0)[0])
            else:
                print "Error: Read files must be named <name>_forward.<filetype> or <name>_reverse.<filetype>"
                exit()

        inputs = zip(forwards_reads, reverse_reads, names)


        threads = 1
        if not args.threads is None and args.threads > 1:
            threads = args.threads

        parallel(runInstance, args, [ProgramRunner("pear", [forwards, reverse, name, threads],
                                                   {"exists": [forwards, reverse]}) for forwards, reverse, name in inputs])
        pool.close()
        pool.join()
        # collect output files
        out_files = []
        out_file_patterns = ["assembled", "discarded", "unassembled.forward", "unassembled.reverse"]
        for forwards, reverse, name in inputs:
            for file_pattern in out_file_patterns:
                out_files.append(glob.glob("%s.%s.%s" % (name, file_pattern, '*'))[0])
        return out_files
        # printVerbose("\t%s sequences assembled, %s contigs discarded, %s sequences discarded" % (-1, -1, -1))
    except KeyboardInterrupt:
        pool.terminate()


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
        parallel(runInstance, args, [ProgramRunner("trimmomatic", [args.forward, args.reverse, args.bdiffs, args.pdiffs,
                                                                   args.oligos, args.processors],
                                                   {"exists": [args.outputFile, args.inputFile, args.oligos],
                                                    "positive": [args.processors],
                                                    "non-Negative": [args.bdiffs, args.pdiffs]})
                                     ])
        return["%s.trim"]
    except KeyboardInterrupt:
        pool.terminate()


def splitOnBarcodes(args, pool=Pool(processes=1)):
    """Splits a fasta/fastq file on a set of barcodes.  An output file will be created for each sample, listing all members
        from that sample.
    :param args: An argparse object with the following parameters:
                    inputFile  File to split
                    barcodes    Tab delimited files of barcodes and their samples
                    outdir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # TODO MAHDI  Comment the rest of the program to replace the pipeline with
    # 1- split_libraries_fastq.py
    # 2- usearch -fastq_filter
    try:
        makeDir(args.outdir)
        printVerbose("Splitting based on barcodes")
        parallel(runInstance, args, [ProgramRunner("barcode.splitter",
                                                   [args.inputFile, args.barcodes,
                                                    os.path.join(args.outdir, "splitOut_")],
                                                   {"exists": [args.inputFile]})
                                     ])
        printVerbose("Demuxed sequences.")

    except KeyboardInterrupt:
        pool.terminate()


def trim(args, pool=Pool(processes=1)):
    """Trims the adapter and barcode from each sequence.

    :param args: An argparse object with the following parameters:
                    inputFasta  Fasta file with sequences to be trimmed
                    oligos      A mothur oligos file with barcodes and primers.  See:
                                    <http://www.mothur.org/wiki/Trim.seqs#oligos>
                    outdir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        printVerbose("Trimming barcodes and adapters")
        makeDir(args.outdir)
        parallel(runInstance, args, [ProgramRunner("trim.seqs",
                                                   [args.inputFasta, args.oligos],
                                                   {"exists": [args.inputFasta, args.oligos]})
                                     ])
        printVerbose("Trimmed sequences.")
        listOfSamples = glob.glob(os.path.join(args.outdir, "splitOut_*"))
    except KeyboardInterrupt:
        pool.terminate()


# TODO decide whether to use this function or splitOnBarcodes.  This one uses SeqIO, the above users fastx
def splitFile(args, pool=Pool(processes=1)):
    """Splits a fastq file on a set of barcodes.  An output file will be created for each sample, listing all members
        from that sample.
    :param args: An argparse object with the following parameters:
                    name        Run Id
                    input_f     Forward Fastq Reads
                    input_r     Reverse Fastq Reads
                    barcodes    Tab delimited files of barcodes and their samples
                    outdir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # Split the cleaned file resulting from joinReads into a user defined
    # number of chunks
    try:
        makeDir(args.outdir)
        splitFileBySample(args.inputFasta, args.groups, args.outdir)
        printVerbose(
            "\tDone splitting file")
        # TODO: eventually send a param to Program running, prevent it from starting after CTRL+C has been invoked
    except KeyboardInterrupt:
        pool.terminate()

def alignSeqs(args, pool=Pool(processes=1)):
    pass



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
    makeDir(args.outdir)
    try:
        if args.program == "macse":
            printVerbose("\t %s Aligning reads using MACSE")
            parallel(runInstance, args, [ProgramRunner("macse_align",
                                                       [args.db, args.db, os.path.join(args.samplesDir, sample)] + [
                                                           os.path.join(args.outdir, sample)] * 3
                                                       , {"exists": []}) for sample in os.listdir(args.samplesDir)])
    except KeyboardInterrupt:
        pool.terminate()


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
    try:
        printVerbose("\t %s Processing MACSE alignments")
        parallel(runInstance, args, [ProgramRunner("macse_format",
                                                   [os.path.join(args.outdir, sample + "_NT"),
                                                    os.path.join(args.outdir, sample + "_AA_macse.fasta"),
                                                    os.path.join(args.outdir, sample + "_NT_macse.fasta"),
                                                    os.path.join(args.outdir, sample + "_macse.csv")],
                                                   {"exists": []}) for sample in os.listdir(args.samplesDir)])

        printVerbose("\tCleaning MACSE alignments")
        # TODO Ask Mahdi what to do with this.  Is this separate step?
        # Remove the reference sequences from the MACSE files and remove the non nucleotide characters from the sequences.
        # we need the datbase seq. names to remove them from the results files

        # TODO:IMPORTANT: Merge the files before doing this.

        dbSeqNames = SeqIO.to_dict(SeqIO.parse(args.db, "fasta")).keys()
        good_seqs = []
        samplesList = os.listdir(args.samplesDir)
        print "Will be processing %s samples " % len(samplesList)
        i = 0
        for sample in samplesList:
            nt_macse_out = os.path.join(args.outdir, sample + "_NT_macse.fasta")
            for mySeq in SeqIO.parse(nt_macse_out, 'fasta'):
                if mySeq.id not in dbSeqNames:
                    mySeq.seq = Seq(str(mySeq.seq[2:]).replace("-", ""))  # remove the !! from the beginning
                    good_seqs.append(mySeq)
            print "completed %s samples" % i
            i += 1
        SeqIO.write(good_seqs, open(os.path.join(args.outdir, "MACSE_OUT_MERGED.fasta"), 'w'), 'fasta')

        printVerbose("\t%s sequences cleaned, %s sequences retained, %s sequences discarded" % (1, 1, 1))

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

            parallel(runInstance, args, [ProgramRunner("chmimera.uchime",
                                                       [args.inputFile, referenceString],
                                                       {"exists": [args.inputFasta, refFile]})
                                         ])
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
    #identify the input file type
    inputFileType = ""
    inputFile = ""
    try:
        # Remove from the accon sequences from the input file
        parallel(runInstance, args, [ProgramRunner("remove.seqs", [args.accnosFile, args.inputFile],
                                                   {"exists": [args.accnosFile, inputFile]})
                                     ])
        # Get the output file name with no directory prefix
        inputFileName = os.path.basename(inputFile)
        splitInputFileName = inputFileName.split(".")
        splitInputFileName.insert(-1, "pick")
        pickOutFile = (".").join(splitInputFileName)

        # TODO move output '*.pick.filetype' file to 'outputDir/*.pick.filetype'
        move("%s/%s" % (os.path.dirname(inputFile), pickOutFile),
             "%s/%s" % (args.outdir,pickOutFile))
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

    #identify the input file type
    # TODO Not sure if mothur actually supports these different imput formats for this command.  Documentation is vague.
    try:
        optionString = mothur_buildOptionString(args, mustFilter=True)
        parallel(runInstance, args, [ProgramRunner("screen.seqs", [args.inputfile, optionString],
                                                   {"exists": [args.inputfile]})
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
        parallel(runInstance, args, [ProgramRunner("make.fastq", [args.inputFasta, args.inputQual],
                                                   {"exists": [args.inputFasta, args.inputQual]})
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
        parallel(runInstance, args, [ProgramRunner("make.fasta", [args.inputFastq], {"exists": [args.inputFastq]})
                                     ])
    except KeyboardInterrupt:
        pool.terminate()


def trimmomatic(args, pool=Pool(processes=1)):
    """Uses a sliding window to identify and trim away areas of low quality.

    :param args: An argparse object with the following parameters:
                    phred	    A boolean toggle.  True for phred-33 scoring, False for phred-64 scoring.
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
        qualScoring = "phrad64"
        if args.phread:
            qualScoring = "phred33"

        parallel(runInstance, args, [ProgramRunner("trimmomatic",
                                                   [qualScoring, args.inputFile, args.outputFile, args.windowSize,
                                                    args.quality, args.minLen],
                                                   {"exists": [args.outputFile, args.inputFile],
                                                    "positive": [args.windowSize, args.quality, args.minLen]})
                                     ])
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
        validate({"exists":[args.aln_out_file, args.caln_userout_file]})

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
#========================================================================================================
#========================================================================================================
#========================================================================================================
def dropShort(args,poo=Pool(processes=1)):
    good_seqs = []
    for seq in SeqIO.parse(args.inputFasta, "fasta"):
        if len(seq.seq) >= int(args.minLenght):
            good_seqs.append(seq)
        else:
            print "seq %s too short (%s bases)" % (seq.id, len(seq.seq))


def serialRename(args, pool=Pool(processes=1)):
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
        parallel(runInstance, args, [
            ProgramRunner("fastx_renamer", [args.input_f, rename_outFile_f], {"exists": [args.input_f]}),
            ProgramRunner("fastx_renamer", [args.input_r, rename_outFile_r], {"exists": [args.input_r]}),
        ])
        printVerbose("\tRenamed %s sequences")
    except KeyboardInterrupt:
        pool.terminate()


# Todo remove this
def doNothing(args, pool=Pool(processes=1)):
    if args.program == "1" and args.arg1 == None:
        print "Error, option 1 requires arg1"
    pass
