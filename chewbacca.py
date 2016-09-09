import argparse
from chewbaccaFunctions import *
import signal


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
        pear + renameSerially.py
        mothur make.contigs + renameSerially.py

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
        vsearch derep_full_length

    9.5 chimeras
        mothur uchime.chimeras

    10 biocode (closed ref)
        vsearch usearch_global

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
    parser.add_argument('--debugtest',default=False)
    subparsers = parser.add_subparsers(dest='action', help='Available commands')


    # ===================================
    # ==  1 Assemble Reads using pear  ==
    # ===================================
    # "pear": programPaths["PEAR"] + " -f \"%s\" -r \"%s\" -o \"%s\" -j %d "
    parser_assemble = subparsers.add_parser('assemble', description="Given a pair of left and right fasta/fastq reads, \
                            or a pair of folder containing the left and right fasta/fastq files, assembles the left \
                            and right read files into contiguous sequence files.  Forwards reads filenames should \
                            end in '_forward.<ext>' or '_R1.<ext>'.  Reverse reads filenames should end in \
                            '_reverse.<ext>' or '_R2.<ext>'.  (where <ext> is the file extension).")
    parser_assemble.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads file or folder.")
    parser_assemble.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads file or folder.")
    parser_assemble.add_argument('-n', '--name',    required=True, help="Assembled File Prefix.")
    parser_assemble.add_argument('-o', '--outdir',  required=True, help="Directory where outputs will be saved.")
    parser_assemble.add_argument('-p', '--pearthreads', type=int, required=False, default=1, help="The number of threads \
                            to use per pear process (default is 1")
    parser_assemble.set_defaults(func=assemble_pear)


    # ======================================
    # ==  2  Split by barcode with FastX  ==
    # ======================================
    # "barcode.splitter": "cat \"%s\" | " + programPaths["FASTX"] + "fastx_barcode_splitter.pl  --bcfile \"%s\" \
    #                                    -prefix \"%s\" --suffix .fastq --bol --mismatches 1",
    parser_demux = subparsers.add_parser('demux_samples', description ="Given a single barcodes file, and a \
                            fasta/fastq file or a folder containing fasta/fastq files, splits each input fasta/fastq \
                            file into separate sample files (based on each sequence's barcode), for each barcode in \
                            the barcodes file.  Sample files of size zero are ignored.")
    parser_demux.add_argument('-i', '--input', required=True, help="Input fasta/fastq file or folder.")
    parser_demux.add_argument('-b', '--barcodes', required=True,
                              help="Tab delimted files of barcodes and corresponding samples.")
    parser_demux.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_demux.set_defaults(func=splitOnBarcodes)


    # ====================================================
    # ==  3 Rename reads serially with renameSequences  ==
    # ====================================================
    # renameSequences(input, output)
    parser_rename = subparsers.add_parser('rename', description="Given a fasta/fastq file or directory of fasta/fastq \
                            files, serially renames each sequence in each file serially, with the filename as a prefix.")
    parser_rename.add_argument('-i', '--input', required=True, help="Input fasta/fastq file or folder.")
    parser_rename.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_rename.add_argument('-f', '--filetype', required=True, help="The filetype of the input files.  Either \
                            'fasta' or 'fastq'.")
    parser_rename.add_argument('-c', '--clip', required=False, default=True, help="True if input file names \
                                contain trailing demux_seqs identifiers.  e.g. True if file name contains '_0', '_1', \
                                '_2', etc..")
    parser_rename.set_defaults(func=renameSequences)


    # ===================================================
    # ==  4 Trims barcodes and adapters using flexbar  ==
    # ===================================================
    # "flexbar":  "flexbar -r \"%s\" -t \"%s\" -ae \"%s\" -a \"%s\""
    parser_trim = subparsers.add_parser('trim_adapters', description="Given a single barcodes file, a single adapters \
                            file, and a single fasta/fastq file or a folder containing fasta/fastq files, removes the \
                            specified adapters (and preceeding barcodes) from all sequences in the given fasta/fastq \
                            files.")
    parser_trim.add_argument('-i', '--input', required=True, help="Input fasta/fastq file or folder.")
    parser_trim.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_trim.add_argument('-a', '--adapters', required=True, help="Forwards Adapters file.")
    parser_trim.add_argument('-arc', '--adaptersrc',  required=True, help="Reverse Complimented Adapters file.")
    parser_trim.add_argument('-u', '--allowedns', required=False, default=0, type=int, help="The number of unknown 'N' \
                            bases a sequence is allowed before being thrown out.  Default: 0.")
    parser_trim.set_defaults(func=trim_flexbar)



    # ======================================
    # ==  5 Clean Reads with Trimmomatic  ==
    # ======================================
    # Clean low-quality reads with trimmomatic
    # "trimomatic":       "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
    # -phred33 input output_cleaned.fastq SLIDINGWINDOW:%windowsize:%minAvgQuality MINLEN:%minLen"
    parser_trimmomatic = subparsers.add_parser('clean_seqs', description="Given a single fastq file or a folder \
                            containing fastq files, trims away areas of low read quality (specified by -q) in each \
                            sequence.  Then, the longest remaining segment (of at least length -m) is recorded as the \
                            new sequence.")
    parser_trimmomatic.add_argument('-i', '--input', required=True, help="Input fastq file or folder")
    parser_trimmomatic.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_trimmomatic.add_argument('-m', '--minLen', type=int, default=200,
                                            help="Minimum length for cleaned sequences")
    parser_trimmomatic.add_argument('-w', '--windowSize', type = int, default=5,
                                            help="Size of the sliding window")
    parser_trimmomatic.add_argument('-q', '--quality', type = int, default=25,
                                            help="Minimum average quality for items in the sliding window")
    parser_trimmomatic.set_defaults(func=trimmomatic)


    # ============================================
    # ==  6 Dereplicate sequences with usearch  ==
    # ============================================
    # "usearch": programPaths["USEARCH"] + " -derep_fulllength \"%s\" -output \"%s\" -uc \"%s\"",
    parser_derep = subparsers.add_parser('dereplicate_fasta', description="Given an fasta file or folder, removes \
                            identical duplicates and subsequences WITHIN EACH FILE, keeping the longest sequence, and \
                            renames it with the number of duplicates as '<longest_sequence_name>_<duplicate count>'.")
    parser_derep.add_argument('-i', '--input', required=True, help="Input fasta file or folder of fasta files.")
    parser_derep.add_argument('-o', '--outdir',  required=True, help="Directory where outputs will be saved.")
    parser_derep.set_defaults(func=u_dereplicate)

    # ============================================
    # ==  6 Dereplicate sequences with vsearch  ==
    # ============================================
    #" vsearch --threads %d --derep_fulllength %s --sizeout --fasta_width 0 --output %s -uc %s",
    parser_vderep = subparsers.add_parser('dereplicate_vsearch', description="Given an fasta file or folder, removes \
                            identical duplicates and subsequences WITHIN EACH FILE, keeping the longest sequence, and \
                            renames it with the number of duplicates as '<longest_sequence_name>_<duplicate count>'.")
    parser_vderep.add_argument('-i', '--input', required=True, help="Input fasta file or folder of fasta files.")
    parser_vderep.add_argument('-o', '--outdir',  required=True, help="Directory where outputs will be saved.")
    parser_vderep.add_argument('-p', '--threads', required=False, type=int, default=2,
                            help="Number of threads to use per query process.")
    parser_vderep.set_defaults(func=dereplicate)



    # ==============================================
    # ==  7 Partition fastas with splitKperFasta  ==
    # ==============================================
    # splitK(inputFasta, prefix, nbSeqsPerFile, filetype):
    parser_split = subparsers.add_parser('partition', description="Given a fasta/fastq file, or folder containing \
                            fasta/fastq files, splits each file into a set of sub files labeled \
                            <filename>_part_x.<ext>, where each subfile contains at most <--chunksize> sequences. \
                            (Where <ext> is the file extension of the original file.")
    parser_split.add_argument('-i', '--input', required=True, help="Input fasta/fastq file or folder.")
    parser_split.add_argument('-o', '--outdir',  required=True, help="Directory where outputs will be saved.")
    parser_split.add_argument('-c', '--chunksize', type=int, required=True, help="Chunksize. (Max number of sequences \
                            per file).")
    parser_split.add_argument('-f', '--filetype', required=True, help="The filetype of the input files.  Either \
                            'fasta' or 'fastq'.")
    parser_split.set_defaults(func=partition)

    # =======================
    # ==  9 Cat the files  ==
    # =======================
    # joinFiles(input_file_list, output_file)
    parser_cat = subparsers.add_parser('merge_files', description="Given a folder containing any type of file, \
                            concatenates all files in that folder together, regardless of their filetype.  Note: Any \
                            file headers/footers will be written as well.")
    parser_cat.add_argument('-i', '--input', required=True, help="Input directory containing files to concatenate.")
    parser_cat.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_cat.add_argument('-n', '--name', required=True, help="Prefix name for the merged file.")
    parser_cat.add_argument('-f', '--fileext', required=True, help="File extension for the output file.")
    parser_cat.set_defaults(func=merge)

    # ===========================
    # ==  9.5 ungap the files  ==
    # ===========================
    # ungap(file_to_clean, output_file_name, gap_char, file_type)
    parser_cat = subparsers.add_parser('ungap_fasta', description="Given a fasta/fastq file or a folder containing \
                            fasta/fastq files, removes alignment characters (specified by -g).")
    parser_cat.add_argument('-i', '--input', required=True, help="Input directory containing files to concatenate.")
    parser_cat.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_cat.add_argument('-f', '--fileext', required=True, help="File extension for the output file.  Either \
                            'fasta', or 'fastq'.")
    parser_cat.add_argument('-g', '--gapchar', required=True, help="A string of one or more characters to remove from \
                            the sequences (but not sequence names) in the input files.")
    parser_cat.set_defaults(func=ungapFasta)

    # ==========================================
    # ==  10 Cluster using vsearch and swarm  ==
    # ==========================================
    # "vsearch": program_paths["VSEARCH"] + "--derep_fulllength \"%s\" --sizeout --fasta_width 0  \
    #
    parser_align = subparsers.add_parser('cluster_seqs', description="Given a fasta file or a folder containing fasta \
                            files, performs clustering on each input file individually.  Outputs a .names file listing \
                            unique representative sequences from each cluster, and the names of the sequences they \
                            represent.  Also outputs a <input_file_name>_seeds.fasta file of the unique representatives.")
    parser_align.add_argument('-i', '--input', required=True, help="Input fasta file/folder")
    parser_align.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_align.add_argument('-n', '--namesfile', required=False, help="A .names file to update.")
    parser_align.add_argument('-s', '--stripcounts', required=False, type=bool, default=True, help="If True, strip \
                            counts from sequence names before clustering.  This allows for the recognition of \
                            sequence names.")
    parser_align.set_defaults(func=cluster)


    # ===================================================
    # ==  12 Closed Ref Picking with fasta references  ==
    # ===================================================
    # --usearch_global  ../9_p_uchime/seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
    #	--userfields query+target+id+alnlen+qcov --userout out  --alnout alnout.txt
    parser_biocode = subparsers.add_parser('query_fasta', description="Given a fasta file, or folder containing \
                            fasta files, aligns sequences asgainst the BIOCODE database.")
    parser_biocode.add_argument('-i', '--input', required=True, help="Input file/folder with fasta files")
    parser_biocode.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_biocode.add_argument('-r', '--referencefasta', required=True, help="Filepath to the curated fasta file to \
                            use as a reference.")
    parser_biocode.add_argument('-x', '--taxinfo', required=True, help="Filepath to a two-column, tab-delimited file \
                            mapping a sequence's fasta id (in the referencefasta file) to a taxonomic identification.")
    parser_biocode.add_argument('-s', '--simmilarity', required=False, default=97, type=int, help="Minimum %  \
                            simmilarity (integer between 0 and 100) between query and reference sequences required for \
                            positive identification. Default: 97")
    parser_biocode.add_argument('-c', '--coverage', required=False, default=85, type=int, help="Minimum % coverage \
                            (integer between 0 and 100) required query and reference sequences required for positive \
                            identification. Default: 85")
    parser_biocode.add_argument('-p', '--threads', required=False, type=int, default=2,
                            help="Number of threads to use per query process.")
    parser_biocode.set_defaults(func=query_fasta)


    # =======================================
    # ==  13 Closed Ref Picking with NCBI  ==
    # =======================================
    # --usearch_global ../9_p_uchime/seeds.pick.fasta  --db /home/mahdi/refs/COI_DOWNLOADED/COI.fasta -id 0.9 \
    #          --userfields query+target+id+alnlen+qcov --userout out --alnout alnout.txt --userfields query+target+id+alnlen+qcov
    parcer_ncbi = subparsers.add_parser('query_ncbi')
    parcer_ncbi.add_argument('-i', '--input', required=True, help="Input Fasta File to clean")
    parcer_ncbi.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parcer_ncbi.add_argument('-p', '--threads', required=False, type=int, default=2,
                             help="Number of threads to use per query process.")
    parcer_ncbi.set_defaults(func=queryNCBI)


    # ========================
    # ==  xx Build Matrix  ==
    # ========================
    parser_build_matrix = subparsers.add_parser('build_matrix', description="Given a single barcodes file with all possible \
                                sample names, a list of the latest names file(s), and a list of initial groups files \
                                (mapping each original, undereplicated sequence to its sample name), builds an OTU \
                                table..  ")
    parser_build_matrix.add_argument('-g', '--groups', required=True, help="Input groups file or folder.")
    parser_build_matrix.add_argument('-n', '--names', required=True, help="Input names file or folder.")
    parser_build_matrix.add_argument('-b', '--barcodes', required=True, help="Input barcodes file.")
    parser_build_matrix.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_build_matrix.set_defaults(func=build_matrix)


    # ==========================
    # ==  xx Annotate Matrix  ==
    # ==========================
    parser_annotate_matrix = subparsers.add_parser('annotate_matrix', description="Given a tabular file mapping sequence IDs \
                                to taxonomic names, and an OTU matrix, renames the identifiable sequence IDs with \
                                taxonomic names.")
    parser_annotate_matrix.add_argument('-i', '--input', required=True, help="Input matrix file or folder of matrix files.")
    parser_annotate_matrix.add_argument('-a', '--annotation', required=True, help="File mapping sequence IDs to taxonomic names.")
    parser_annotate_matrix.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_annotate_matrix.set_defaults(func=annotate_matrix)


    # =================================
    # ==  xx Convert fastq to fasta  ==
    # =================================
    # translateFastqToFasta(inputFastQ, outputFasta):
    parser_toFasta = subparsers.add_parser('convert_fastq_to_fasta')
    parser_toFasta.add_argument('-i', '--input', required=True, help="Fastq file or folder containing fastq files to "
                                 "translate")
    parser_toFasta.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_toFasta.set_defaults(func=makeFasta)

    global args
    args, unknown = parser.parse_known_args()
    if unknown:
        print "\nIgnoring unknown args: " + ', '.join(['%s']*len(unknown))% tuple(unknown)
    if args.verbose:
        logging.basicConfig(format=FORMAT, level=logging.DEBUG, datefmt=DATEFMT)
    else:
        logging.basicConfig(format=FORMAT, level=logging.ERROR, datefmt=DATEFMT)

    printVerbose.VERBOSE = args.verbose
    logging.debug("Initial ARGS are: %s", args)
    print("\t\t")
    dryRun = args.dryRun
    signal.signal(signal.SIGTSTP, signal.SIG_IGN)
    args.func(args, args.debugtest)


if __name__ == "__main__":
    main(sys.argv)
