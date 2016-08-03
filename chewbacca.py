import argparse
from classes.Helpers import *
from chewbaccaFunctions import *
from multiprocessing import Pool



# Program Version
version = "0.01"
# Time format for printing
FORMAT = "%(asctime)s  %(message)s"
# Date format for printing
DATEFMT = "%m/%d %H:%M:%S"

# TODO add output file option to all functions, and rename/move all the output files.
# rm -r rslt/;python chewbacca.py preprocess -n test -f ~/ARMS/testARMS/testData/20K_R2.fq -r ~/ARMS/testARMS/testData/20K_R1.fq -b /home/greg/ARMS/testARMS/testData/barcodes.txt -o rslt

"""
Supported operations:
    1. assemble
        pear + rename_sequences.py
        mothur make.contigs + rename_sequences.py

    2 demux
        fastx barcode splitter
        mothur split.groups

    3 trim (barcodes/adapters)
        mothur trim.seqs
        flexbar

    4 clean (sliding window)
        mothur trim.seqs
        trimmomatic

    5 dereplicate
        mothur unique.seqs
        ~/ARMS/bin/usearch7.0.1090

    6 align
        mothur align.seqs

    7 split
        splitKperFasta.py

    8 macsealign
        macse

    9 cluster
        macse
        ~/bin/vsearch/bin/vsearch

    9.5 chimeras
        mothur uchime.chimeras


    10 biocode (closed ref)
        ~/bin/vsearch/bin/vsearch-1.1.1-linux-x86_64

    11 NCBI    (closed ref)
        ~/bin/vsearch/bin/vsearch-1.1.1-linux-x86_64

    12 buildmatrix
        mothur Make.biom
        ~/bin/builMatrix.py
"""

# TODO parser args should probably live next to their program strings (I.E. in Program Runner or another class)
def main(argv):
    """Parses command line args, builds an argparse.ArgumentParser, and runs the chosen command.
        Otherwise, prints usage.

    :param argv: Command line arguments as a list of strings
    """
    parser = argparse.ArgumentParser(description="arms description", epilog="arms long description")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    parser.add_argument("--verbose", default=True, help="increase output verbosity", action="store_true")
    parser.add_argument('-t', '--threads', type=int, default=1)
    parser.add_argument('--dryRun', default=False)
    subparsers = parser.add_subparsers(dest='action', help='Available commands')

    # ===========================================
    # ==  Assemble Reads using mothur or pear  ==
    # ===========================================
    parser_assemble = subparsers.add_parser('assemble')
    parser_assemble.add_argument('-p', '--program', required=True, help="The number of threads to use")
    parser_assemble.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads")
    parser_assemble.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads")
    parser_assemble.add_argument('-n', '--name',    required=True, help="Assembled File Prefix")
    parser_assemble.add_argument('-o', '--outdir',  required=True, help="Directory where outputs will be saved")
    parser_assemble.set_defaults(func=assemble)

    # for mothur assembler
    mothur_assembler = parser_assemble.add_argument_group('mothur', 'Mothur options')
    mothur_assembler.add_argument('-g', '--oligos',  help="Oligos file with barcode and primer sequences")
    mothur_assembler.add_argument('-bd', '--bdiffs', help="# of allowed barcode mismatches")
    mothur_assembler.add_argument('-pd', '--pdiffs', help="# of allowed primer mismatches")

    # for pear assembler
    pear_assembler = parser_assemble.add_argument_group('pear', 'Pear options')
    pear_assembler.add_argument('-m', '--maxlength', type=int, help="Maximum length for assembled sequences")
    pear_assembler.add_argument('-t', '--threads',   type=int, help="The number of threads to use")




    '''
    # Assemble reads
    # "pear":             programPaths["PEAR"] + " -f \"%s\" -r \"%s\" -o \"%s\" -j \"%s\" -m %d",
    # "make.contigs": "mothur \'#make.contigs(ffastq=%s, rfastq=%s, bdiffs=1, pdiffs=2, oligos=%s, processors=%s)\'"
    parser_assemble = subparsers.add_parser('assemble')
    parser_assemble.add_argument('-p', '--program', type=int, required=True, help="The number of threads to use")
    parser_assemble.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads")
    parser_assemble.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads")
    parser_assemble.add_argument('-n', '--name', required=True, help="Assembled File Prefix")
    parser_assemble.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_assemble.add_argument('-t', '--threads', type=int, help="The number of threads to use")
    # options
    # for mothur
    parser_assemble.add_argument('-g', '--oligos', required=True,
                                    help="Oligos file with barcode and primer sequences")
    parser_assemble.add_argument('-bd', '--bdiffs', required=True, help="# of allowed barcode mismatches")
    parser_assemble.add_argument('-pd', '--pdiffs', required=True, help="# of allowed primer mismatches")
    # for pear
    parser_assemble.add_argument('-m', '--maxLen', type=int,
                                 help="Maximum length for assembled sequences")
    parser_assemble.set_defaults(func=assemble)
    '''

    # Rename reads serially with Fastx
    # "fastx_renamer":    programPaths["FASTX"] + "fastx_renamer -n COUNT -i \"%s\" -o \"%s\" -Q 33",
    parser_serialize = subparsers.add_parser('rename')
    parser_serialize.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads")
    parser_serialize.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads")
    parser_serialize.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_serialize.set_defaults(func=serialRename)



    # Trims the barcode and adapters/primers with Mothur
    # "trim.seqs":        "mothur \'#trim.seqs(fasta=\"%s\", oligos=\"%s\", maxambig=1, maxhomop=8, \
    #                                minlength=300, maxlength=550, bdiffs=1, pdiffs=2)\'",
    parser_mothurTrim = subparsers.add_parser('mothurTrim')
    parser_mothurTrim.add_argument('-i', '--inputFasta', required=True, help="Input Fasta/Fastq File")
    parser_mothurTrim.add_argument('-g', '--oligos', required=True, help="Mothur oligos file")
    parser_mothurTrim.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_mothurTrim.set_defaults(func=trim)



    """
    # Assemble reads with Pear
    # "pear":             programPaths["PEAR"] + " -f \"%s\" -r \"%s\" -o \"%s\" -j \"%s\" -m %d",
    # "make.contigs": "mothur \'#make.contigs(ffastq=%s, rfastq=%s, bdiffs=1, pdiffs=2, oligos=%s, processors=%s)\'"
    parser_assemble = subparsers.add_parser('assemble')
    parser_assemble.add_argument('-n', '--name', required=True, help="Run Id")
    parser_assemble.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads")
    parser_assemble.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads")
    parser_assemble.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_assemble.add_argument('-m', '--maxLen', type=int,
                                    help="Maximum length for assembled sequences")
    parser_assemble.add_argument('-t', '--threads', type=int, help="The number of threads to use")
    parser_assemble.set_defaults(func=assembleReads)


    # Assembles forward and reverse fastq reads
    # "make.contigs": "mothur \'#make.contigs(ffastq=%s, rfastq=%s, bdiffs=1, pdiffs=2, oligos=%s, processors=%s)\'"
    parser_makeContigs = subparsers.add_parser('makeContigs')
    parser_makeContigs.add_argument('-f', '--forward', required=True, help="Forward read fastq file")
    parser_makeContigs.add_argument('-r', '--reverse', required=False, help="Reverse read fastq file")
    parser_makeContigs.add_argument('-g', '--oligos', required=True,
                                    help="Oligos file with barcode and primer sequences")
    parser_makeContigs.add_argument('-bd', '--bdiffs', required=True, help="# of allowed barcode mismatches")
    parser_makeContigs.add_argument('-pd', '--pdiffs', required=True, help="# of allowed primer mismatches")
    parser_makeContigs.add_argument('-p', '--procs', required=True, help="Number of processors to use")
    """


    # Split file by barcode with Fastx
    # "barcode.splitter": "cat \"%s\" | " + programPaths["FASTX"] + "fastx_barcode_splitter.pl  --bcfile \"%s\" \
    #                                    -prefix \"%s\" --suffix .fastq --bol --mismatches 1",
    parser_demux = subparsers.add_parser('demux')
    parser_demux.add_argument('-i', '--inputFile', required=True, help="Input fasta/fastq")
    parser_demux.add_argument('-b', '--barcodes', required=True,
                              help="Tab delimted files of barcodes and their samples")
    parser_demux.add_argument('-o', '--outdir', required=True, help="Output directory")
    parser_demux.set_defaults(func=splitOnBarcodes)



    # Split file by samples
    parser_split = subparsers.add_parser('partition')
    parser_split.add_argument('-n', '--name', required=True, help="Run Id")
    parser_split.add_argument('-i', '--inputFasta', required=True, help="Input fasta file to split")
    parser_split.add_argument('-g', '--groups', required=True,
                              help=" Groups file in same format as that generated by mothur")
    parser_split.add_argument('-o', '--outdir', required=True, help="Output directory")
    parser_split.set_defaults(func=splitFile)



    # Align Reads with MACSE
    # "macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
    #                                \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
    #                                -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
    #                                -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",
    parser_align = subparsers.add_parser('align')
    parser_align.add_argument('-s', '--samplesDir', required=True,
                              help="Directory containing the samples file required for clustering")
    parser_align.add_argument('-p', '--program', required=True,
                              help="Program name for cleaning. Available options are: ... ")
    parser_align.add_argument('-d', '--db', required=True, help="Database against which to align and filter reads")
    parser_align.add_argument('-o', '--outdir', required=True, help=" Output directory")
    parser_align.set_defaults(func=macseAlignSeqs)



    # Clean Aligned Reads with MACSE
    # "macse_format":     "java -jar " + programPaths["MACSE"] + "  -prog exportAlignment -align \"%s\" \
    #                                -charForRemainingFS - -gc_def 5 -out_AA \"%s\" -out_NT \"%s\" -statFile \"%s\"",
    parser_macseClean = subparsers.add_parser('macseClean')
    parser_macseClean.add_argument('-s', '--samplesDir', required=True,
                                            help="Directory containing the samples file required for clustering")
    parser_macseClean.add_argument('-o', '--outdir', required=True, help=" Output directory")
    parser_macseClean.set_defaults(func=macseCleanAlignments)



    # Clean low-quality reads with trimmomatic
    # "trimomatic":       "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
    # -phred33 input output_cleaned.fastq SLIDINGWINDOW:%windowsize:%minAvgQuality MINLEN:%minLen"
    parser_trimmomatic = subparsers.add_parser('trimmomatic')
    parser_trimmomatic.add_argument('-i', '--inputFile', required=True, help="Run Id")
    parser_trimmomatic.add_argument('-o', '--outputFile', required=True, help="Forward Fastq Reads")
    parser_trimmomatic.add_argument('-m', '--minLen', required=True, type=int,
                                            help="Minimum length for cleaned sequences")
    parser_trimmomatic.add_argument('-p', '--phred33', required=True, type=bool,
                                            help="T for phred33, F for phred64 quality scores")
    parser_trimmomatic.add_argument('-w', '--windowSize', required=True, help="Size of the sliding window")
    parser_trimmomatic.add_argument('-q', '--quality', required=True,
                                            help="Minimum average quality for items in the sliding window")
    parser_trimmomatic.set_defaults(func=trimmomatic)




    # Find Chimeras with Mothur
    # "chmimera.uchime":  "mothur \'#chimera.uchime(fasta=\"%s\", name=\"%s\")\'",
    parser_chimera = subparsers.add_parser('findChimeras')
    parser_chimera.add_argument('-i', '--inputFasta', required=True, help="Input Fasta File to clean")
    parser_chimera.add_argument('-o', '--outdir', required=True, help="Output directory")
    parser_chimera.add_argument('-p', '--program', required=False, default="uchime",
                                help="Program for detecting and removing chimeras. Default is uchime")
    refType = parser_chimera.add_mutually_exclusive_group(required=True)
    refType.add_argument('-n', '--namesFile', help="Names file to reference")
    refType.add_argument('-r', '--refDB', help="Database file to reference")
    parser_chimera.set_defaults(func=findChimeras)



    # Screen seqs with Mothur
    # outputs: *.good.fasta and *.bad.accnos
    # "screen.seqs": "mothur \'#screen.seqs(fasta=%s, %s)\'",
    parser_screen= subparsers.add_parser('screen')
    parser_screen.add_argument('-i', '--inputfile', required=True, action ='store',
                               help="File path to the file to clean")
    parser_screen.add_argument('-o', '--outdir', help="File path to your output directory")
    # optional filters
    parser_screen.add_argument('--start', help="Maximum allowable sequence starting index")
    parser_screen.add_argument('--end', help="Minimum allowable sequence ending index")
    parser_screen.add_argument('--minlength', help="Minimum allowable sequence length")
    parser_screen.add_argument('--maxlength', help="Maximum allowable sequence length")
    parser_screen.add_argument('--maxambig', help="Maxmimum number of allowed ambiguities")
    parser_screen.add_argument('--maxn', help="Maximum number of allowed N's")
    parser_screen.add_argument('--maxhomop', help="Maximum allowable homopolymer length")

    # optional aux files to update
    parser_screen.add_argument('-g', '--groups', help="Groups file to update")
    parser_screen.add_argument('-n', '--names', help="Names file to update")
    parser_screen.add_argument('-r', '--alnReport', help="Alignment report to update")
    parser_screen.add_argument('-c', '--contigsReport', help="Contigs report to update")
    parser_screen.add_argument('-s', '--summaryFile', help="SummaryFile to update")
    parser_screen.set_defaults(func=screenSeqs)



    # Remove sequences with Mothur
    # "remove.seqs": "mothur \'#remove.seqs(accnos=%s, %s)\'",
    parser_chimera = subparsers.add_parser('removeSeqs')
    parser_chimera.add_argument('-a', '--accnosFile', required=True, help="Accnos File listing sequences to remove")
    parser_chimera.add_argument('-o', '--outdir', required=True, help="Output directory")
    refType = parser_chimera.add_mutually_exclusive_group(required=True)
    refType.add_argument('-f', '--fasta', help="Fasta/Fastq file to clean")
    refType.add_argument('-l', '--list', help="List file to clean")
    refType.add_argument('-g', '--groups', help="Groups file to clean")
    refType.add_argument('-n', '--names', help="Names file to clean")
    refType.add_argument('-c', '--count', help="Count file to clean")
    refType.add_argument('-r', '--alnReport', help="Alignment report to clean")
    parser_chimera.set_defaults(func=removeSeqs)


    # Drop short reads
    parser_dropShort = subparsers.add_parser('dropShort')
    parser_dropShort.add_argument('-n', '--name', required=True, help="Run Id")
    parser_dropShort.add_argument('-i', '--inputFasta', required=True, help="Clean inputs File")
    parser_dropShort.add_argument('-f', '--namesFile', required=True, help="Updated names file")
    parser_dropShort.add_argument('-l', '--minLenght', required=True, help="Min. length to keep")
    # todo -o is resereved for outdir
    # parser_dropShort.add_argument('-o', '--outFasta', required=True, help="Output file filtered on lenght")
    parser_dropShort.add_argument('-s', '--outNames', required=True, help="Updated names file")
    parser_dropShort.set_defaults(func=dropShort)



    # Drop short reads
    parser_cluster = subparsers.add_parser('cluster-swarm')
    # parser_cluster.set_defaults(func=clusterReads)

    # Prescreen sequences for frameshifts (which will get fed to MACSE)
    # screen(aln, caln)
    parser_chimera = subparsers.add_parser('prescreen')
    parser_chimera.add_argument('-a', '--aln_out_file', required=True, help=".aln output file from vsearch")
    parser_chimera.add_argument('-c', '--caln_userout_file', required=True, help="caln output file from vsearch.")
    parser_chimera.add_argument('-o', '--outdir',  required=True, help="Directory where outputs will be saved")
    parser_chimera.set_defaults(func=prescreen)

    # Convert fastq to fasta
    # "make.fasta":       "mothur \'#fastq.info(fastq=%s,fasta=T)\'"
    parser_toFasta = subparsers.add_parser('makeFasta')
    parser_toFasta.add_argument('-i', '--inputFastq', required=True, help="Input Fastq File")
    parser_toFasta.set_defaults(func=makeFasta)



    # Convert fasta to fastq
    # "make.fastq":       "mothur \'#make.fastq(fasta=%s,qfile=%s)\'",
    parser_toFastq = subparsers.add_parser('makeFastq')
    parser_toFastq.add_argument('-i', '--inputFasta', required=True, help="Input Fasta File")
    parser_toFastq.add_argument('-q', '--inputQual', required=True, help="Input qual File")
    parser_toFastq.set_defaults(func=makeFastq)

    # Todo remove this
    testParser = subparsers.add_parser('test')
    testParser.add_argument('-p', '--program',  help="Input Fasta File")
    group1 = testParser.add_argument_group("group1","G1help")
    group1.add_argument('-x', '--arg1', help="arg1 text")
    group2 = testParser.add_argument_group("group2", "G1help")
    group2.add_argument('-y', '--arg2', help="arg2 text")
    testParser.set_defaults(func=doNothing)

    global args, pool
    args, unknown = parser.parse_known_args()
    if unknown:
        print "\nIgnoring unknown args: " + ', '.join(['%s']*len(unknown))% tuple(unknown)
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
