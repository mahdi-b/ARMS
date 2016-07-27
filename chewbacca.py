import argparse
import glob
import logging
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from classes.Helpers import printVerbose
from classes.Helpers import makeDirOrdie
from classes.Helpers import splitFileBySample
from classes.ProgramRunner import ProgramRunner
from multiprocessing import Pool
from shutil import copyfile

# Program Version
version = "0.01"
# Time format for printing
FORMAT = "%(asctime)s  %(message)s"
# Date format for printing
DATEFMT = "%m/%d %H:%M:%S"

# rm -r rslt/;python chewbacca.py preprocess -n test -f ~/ARMS/testARMS/testData/20K_R2.fq -r ~/ARMS/testARMS/testData/20K_R1.fq -b /home/greg/ARMS/testARMS/testData/barcodes.txt -o rslt

def runInstance(args, myInstance):
    """Runs an instance of a ProgramRunner.  Calls ProgramRunner.run() for a ProgramRunner object.'
        :param myInstance A fully initalized ProgramRunner object to run.
    """
    # Use the global version to facilitate calling workers
    print "runinstance"
    if args.dryRun:
        logging.info(myInstance.dryRun())
    else:
        #logging.info(myInstance.dryRun())
        myInstance.run()

def parallel(function, args, runners, pool=Pool(processes=1)):
    '''
    Executes one or more ProgramRunners in parallel.
    :param function:    The function to call over ProgramRunners.  Generally runInstance().
    :param args:        Arguments from the argparse object.
    :param runners:     A list of fully initalized ProgramRunners to execute.
    :param pool:        An initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    :return:
    '''

    data = [(args, runner) for runner in runners]
    for arg, runner in data:
            runInstance(arg,runner)
    #pool.map(function, data)


def serialRename(args, pool=Pool(processes=1)):
    """Renames sequences in a fastq file as 1,2,3,...
    :param args: An argparse object with the following parameters:
                    input_f     Forward Fastq Reads
                    input_r     Reverse Fastq Reads
                    outDir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        # Make the output directory, or abort if it already exists
        makeDirOrdie(args.outDir)
        printVerbose("\tRenaming sequences")
        # "~/programs/fastx/bin/fastx_renamer -n COUNT -i %s %s"
        rename_outFile_f = os.path.join(args.outDir, os.path.basename(args.input_f) + "_renamed")
        rename_outFile_r = os.path.join(args.outDir, os.path.basename(args.input_r) + "_renamed")
        parallel(runInstance, args, [
                                ProgramRunner("fastx_renamer",[args.input_f, rename_outFile_f], {"exists":[args.input_f]}),
                                ProgramRunner("fastx_renamer",[args.input_r, rename_outFile_r], {"exists":[args.input_r]}),
                                ])
        printVerbose("\tRenamed %s sequences")
    except KeyboardInterrupt:
        pool.terminate()

def assembleReads(args, pool=Pool(processes=1)):
    """Assembles reads from two (left and right) fastq files.
    :param args: An argparse object with the following parameters:
                    name        Textual ID for the data set
                    input_f     Forward Fastq Reads
                    input_r     Reverse Fastq Reads
                    threads     The number of threads to use durring assembly.
                    outDir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """

    # *****************************************************************************************
    # Making the contigs using Pear
    # "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s"
    try:
        printVerbose("\tAssembling reads")
        # TODO figure out a better way to handle filenames
        rename_outFile_f = os.path.join(args.outDir, os.path.basename(args.input_f) + "_renamed")
        rename_outFile_r = os.path.join(args.outDir, os.path.basename(args.input_r) + "_renamed")

        assembledPrefix = os.path.join(args.outDir, args.name)
        parallel(runInstance, args, [ProgramRunner("pear",
                                             [rename_outFile_f, rename_outFile_r, assembledPrefix, args.threads],
                                             {"exists":[rename_outFile_f, rename_outFile_r]})
                                ])

        assembledFastqFile = os.path.join(args.outDir, args.name + ".assembled.fastq")
        printVerbose("\t%s sequences assembled, %s contigs discarded, %s sequences discarded" % (-1, -1, -1))
    except KeyboardInterrupt:
        pool.terminate()


def splitOnBarcodes(args, pool=Pool(processes=1)):
    """Splits a fasta/fastq file on a set of barcodes.  An output file will be created for each sample, listing all members
        from that sample.
    :param args: An argparse object with the following parameters:
                    name        Run Id
                    barcodes    Tab delimited files of barcodes and their samples
                    outDir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # TODO MAHDI  Comment the rest of the program to replace the pipeline with
    # 1- split_libraries_fastq.py
    # 2- usearch -fastq_filter
    try:
        # TODO figure out a better way to handle filenames
        assembledFastqFile = os.path.join(args.outDir, args.name + ".assembled.fastq")

        printVerbose("Splitting based on barcodes")
        parallel(runInstance, args, [ProgramRunner("barcode.splitter",
                                             [assembledFastqFile, args.barcodes, os.path.join(args.outDir, "splitOut_")],
                                             {"exists": [assembledFastqFile]})
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
                    outDir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    try:
        makeDirOrdie(args.outDir)
        printVerbose("Splitting based on barcodes")
        parallel(runInstance, args, [ProgramRunner("trim.seqs",
                                                   [args.inputFasta, args.oligos],
                                                   {"exists": [args.inputFasta, args.oligos]})
                                     ])
        printVerbose("Trimmed sequences.")
        listOfSamples = glob.glob(os.path.join(args.outDir, "splitOut_*"))
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
                    outDir      Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # Split the cleaned file resulting from joinReads into a user defined
    # number of chunks
    try:
        makeDirOrdie(args.outDir)
        splitFileBySample(args.inputFasta, args.groups, args.outDir)
        printVerbose(
            "\tDone splitting file")  # TODO: eventually send a param to Program running, prevent it from starting after CTRL+C has been invoked
    except KeyboardInterrupt:
        pool.terminate()


def alignSeqs(args, pool=Pool(processes=1)):
    """Aligns sequences by iteratively adding them to a known good alignment.

     :param args: An argparse object with the following parameters:
                    db                  Database against which to align and filter reads
                    samplesDir          Directory containig the samples to be cleaned
                    outDir              Directory where outputs will be saved
     :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    #"macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
    #                                    \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
    #                                    -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
    #                                    -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",
    makeDirOrdie(args.outDir)
    try:
        if args.program == "macse":
            printVerbose("\t %s Aligning reads using MACSE")
            parallel(runInstance, args, [ProgramRunner("macse_align",
                                                 [args.db, args.db, os.path.join(args.samplesDir, sample)] + [
                                                     os.path.join(args.outDir, sample)] * 3
                                                 , {"exists": []}) for sample in os.listdir(args.samplesDir)])
    except KeyboardInterrupt:
        pool.terminate()


def cleanAlignments(args, pool=Pool(processes=1)):
    """Removes non-nucleotide characters in MACSE aligned sequences for all fasta files in the samples directory
        (the samplesDir argument).
    :param args: An argparse object with the following parameters:
                    samplesDir          Directory containig the samples to be cleaned
                    outDir              Directory where outputs will be saved
    :param pool: A fully initalized multiprocessing.Pool object.  Defaults to a Pool of size 1.
    """
    # "macse_format":     "java -jar " + programPaths["MACSE"] + "  -prog exportAlignment -align \"%s\" \
    #
    #                                  -charForRemainingFS - -gc_def 5 -out_AA \"%s\" -out_NT \"%s\" -statFile \"%s\""
    try:
        printVerbose("\t %s Processing MACSE alignments")
        parallel(runInstance, args, [ProgramRunner("macse_format",
                                             [os.path.join(args.outDir, sample + "_NT"),
                                              os.path.join(args.outDir, sample + "_AA_macse.fasta"),
                                              os.path.join(args.outDir, sample + "_NT_macse.fasta"),
                                              os.path.join(args.outDir, sample + "_macse.csv")],
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
            nt_macse_out = os.path.join(args.outDir, sample + "_NT_macse.fasta")
            for mySeq in SeqIO.parse(nt_macse_out, 'fasta'):
                if mySeq.id not in dbSeqNames:
                    mySeq.seq = Seq(str(mySeq.seq[2:]).replace("-", ""))  # remove the !! from the beginning
                    good_seqs.append(mySeq)
            print "completed %s samples" % i
            i += 1
        SeqIO.write(good_seqs, open(os.path.join(args.outDir, "MACSE_OUT_MERGED.fasta"), 'w'), 'fasta')

        printVerbose("\t%s sequences cleaned, %s sequences retained, %s sequences discarded" % (1, 1, 1))

    except KeyboardInterrupt:
        pool.terminate()


def removeChimeras(args, pool=Pool(processes=1)):
    """

    :param args: A list of arguments to the function
    :param pool: A multiprocessing.Pool object.
    """
    if args.program == "uchime":
        # *****************************************************************************************
        # findind the chimeras
        # "chmimera.uchime": "mothur #chimera.uchime(fasta=%s, name=%s)",
        parallel(runInstance, args, [ProgramRunner("chmimera.uchime",
                                                   [args.inputFile, args.namesFile],
                                                   {"exists": [args.inputFile, args.namesFile]})
                               ])
        # *****************************************************************************************
        # removing the chimeras from input file
        # remove.seqs(accnos=Nucleotidealignment.uchime.accnos, fasta=Nucleotidealignment.fasta)

        uchimeAccnos = glob.glob(os.path.join(os.path.dirname(args.inputFile), "*uchime.chimeras"))[0]
        print "outFile is %s " % uchimeAccnos
        p = ProgramRunner("remove.seqs", [uchimeAccnos, args.inputFile], {"exists": [uchimeAccnos]})
        parallel(runInstance, args, p)

        # *****************************************************************************************
        # Renaming the outfile and updating the names file
        # remove from the names file the sequences that were removed

        pickOutFile = glob.glob(os.path.join(os.path.dirname(args.inputFile), "*pick.fasta"))[0]
        copyfile(pickOutFile, args.outFasta)

        # TO CONITNUE: UPDATE THE NAMES file


    else:
        raise Exception("unknown program %s for chimera detection or removal" % args.program)

    printVerbose("\t removed %s chimeric sequences")


def dropShort(args, pool=Pool(processes=1)):
    """
    :param args: A list of arguments to the function
    :param pool: A multiprocessing.Pool object.
    """
    # Dropping short sequences
    good_seqs = []
    for seq in SeqIO.parse(args.inputFasta, "fasta"):
        if len(seq.seq) >= int(args.minLenght):
            good_seqs.append(seq)
        else:
            print "seq %s too short (%s bases)" % (seq.id, len(seq.seq))


def main(argv):
    """Parses command line args, builds an argparse.ArgumentParser, and runs the chosen command.
        Otherwise, prints usage.

    :param argv: Command line arguments as a list of strings
    """
    parser = argparse.ArgumentParser(description="arms description", epilog="arms long description")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    parser.add_argument("--verbose", default=True, help="increase output verbosity", action="store_true")
    parser.add_argument('-t', '--threads', type=int, default=1)
    parser.add_argument('--dryRun', action='store_true', default=False)
    subparsers = parser.add_subparsers(dest='action', help='Available commands')

    # Rename reads serially
    parser_serialize = subparsers.add_parser('serialize')
    parser_serialize.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads")
    parser_serialize.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads")
    parser_serialize.add_argument('-o', '--outDir', required=True, help="Directory where outputs will be saved")
    parser_serialize.set_defaults(func=serialRename)

    # Trims the barcode and adapters/primers
    parser_trim = subparsers.add_parser('trim')
    parser_trim.add_argument('-i', '--inputFasta', required=True, help="Input Fasta/Fastq File")
    parser_trim.add_argument('-g', '--oligos', required=True, help="Mothur oligos file")
    parser_trim.add_argument('-o', '--outDir', required=True, help="Directory where outputs will be saved")
    parser_trim.set_defaults(func=trim)

    # Assemble reads
    parser_assemble = subparsers.add_parser('assemble')
    parser_assemble.add_argument('-n', '--name', required=True, help="Run Id")
    parser_assemble.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads")
    parser_assemble.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads")
    parser_assemble.add_argument('-o', '--outDir', required=True, help="Directory where outputs will be saved")
    #parser_assemble.add_argument('-t', '--threads', required=True, help="The number of threads to use")
    parser_assemble.set_defaults(func=assembleReads)

    # split file by barcode
    parser_demux = subparsers.add_parser('demux')
    parser_demux.add_argument('-n', '--name', required=True, help="Run Id")
    parser_demux.add_argument('-b', '--barcodes', required=True,
                                   help="Tab delimted files of barcodes and their samples")
    parser_demux.add_argument('-o', '--outDir', required=True, help="Output directory")
    parser_demux.set_defaults(func=splitOnBarcodes)

    # split file by samples
    parser_split = subparsers.add_parser('partition')
    parser_split.add_argument('-n', '--name', required=True, help="Run Id")
    parser_split.add_argument('-i', '--inputFasta', required=True, help="Input fasta file to split")
    parser_split.add_argument('-g', '--groups', required=True,
                              help=" Groups file in same format as that generated by mothur")
    parser_split.add_argument('-o', '--outDir', required=True, help="Output directory")
    parser_split.set_defaults(func=splitFile)

    # Align Reads
    parser_align = subparsers.add_parser('align')
    parser_align.add_argument('-s', '--samplesDir', required=True,
                              help="Directory containing the samples file required for clustering")
    parser_align.add_argument('-p', '--program', required=True,
                              help="Program name for cleaning. Available options are: ... ")
    parser_align.add_argument('-d', '--db', required=True, help="Database against which to align and filter reads")
    parser_align.add_argument('-o', '--outDir', required=True, help=" Output directory")
    parser_align.set_defaults(func=alignSeqs)


    # Clean Aligned Reads
    parser_align = subparsers.add_parser('clean')
    parser_align.add_argument('-s', '--samplesDir', required=True,
                              help="Directory containing the samples file required for clustering")
    parser_align.add_argument('-o', '--outDir', required=True, help=" Output directory")
    parser_align.set_defaults(func=cleanAlignments)


    # Remove Chimeras
    parser_chimera = subparsers.add_parser('removeChimeras')
    parser_chimera.add_argument('-n', '--name', required=True, help="Run Id")
    parser_chimera.add_argument('-i', '--inputFile', required=True, help="Clean inputs File")
    parser_chimera.add_argument('-p', '--program', required=False, default="uchime",
                                help="Program for detecting and removing chimeras. Default is uchime")
    parser_chimera.add_argument('-f', '--namesFile', required=True, help="Updated names file")
    parser_chimera.add_argument('-o', '--outFasta', required=True, help="Updated names file")
    parser_chimera.add_argument('-s', '--outNames', required=True, help="Updated names file")
    parser_chimera.set_defaults(func=removeChimeras)

    # Drop short reads
    parser_dropShort = subparsers.add_parser('dropShort')
    parser_dropShort.add_argument('-n', '--name', required=True, help="Run Id")
    parser_dropShort.add_argument('-i', '--inputFasta', required=True, help="Clean inputs File")
    parser_dropShort.add_argument('-f', '--namesFile', required=True, help="Updated names file")
    parser_dropShort.add_argument('-l', '--minLenght', required=True, help="Min. length to keep")
    parser_dropShort.add_argument('-o', '--outFasta', required=True, help="Output file filtered on lenght")
    parser_dropShort.add_argument('-s', '--outNames', required=True, help="Updated names file")
    parser_dropShort.set_defaults(func=dropShort)

    # Drop short reads
    parser_cluster = subparsers.add_parser('cluster-swarm')
    # parser_cluster.set_defaults(func=clusterReads)


    # Convert fastq to fasta

    parser_toFastq = subparsers.add_parser('toFastq')
    parser_toFastq.add_argument('-i', '--inputFasta', required=True, help="Input Fasta File")
    parser_toFastq.add_argument('-q', '--inputQual', required=True, help=" Input qual file")
    parser_toFastq.set_defaults(func=convert)

    global args, pool
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format=FORMAT, level=logging.DEBUG, datefmt=DATEFMT)
    else:
        logging.basicConfig(format=FORMAT, level=logging.ERROR, datefmt=DATEFMT)

    printVerbose.VERBOSE = args.verbose
    printVerbose("Running with %s process(es)" % args.threads)
    pool = Pool(args.threads)
    logging.debug("Initial ARGS are: %s", args)
    print("\t\t")
    dryRun = args.dryRun
    args.func(args, pool)

    pool.close()
    pool.join()

if __name__ == "__main__":
    main(sys.argv)
