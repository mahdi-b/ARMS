import argparse
import signal

from clean.Clean_Deep_Repair_Command import Clean_Deep_Repair_Command
from assemble.Assemble_Command import Assemble_Command
from classes.Helpers import *
from clean.Clean_Adapters_Command import Clean_Adapters_Command
from clean.Clean_Deep_Command import Clean_Deep_Command
from clean.Clean_Quality_Command import Clean_Quality_Command
from cluster.Cluster_Command import Cluster_Command
from demux.Demux_Command import Demux_Command
from dereplicate.Dereplicate_Command import Dereplicate_Command
from otu.Annotate_OTU_Table_Command import Annotate_OTU_Table_Command
from otu.Build_OTU_Table_Command import Build_OTU_Table_Command
from otu.Query_OTU_DB_Command import Query_OTU_DB_Command
from otu.Query_OTU_Fasta_Command import Query_OTU_Fasta_Command
from preclean.Preclean_Command import Preclean_Command
from rename.Rename_Command import Rename_Command
from util.Convert_Fastq_Fasta_Command import Convert_Fastq_Fasta_Command
from util.Merge_Command import Merge_Command
from util.Partition_Command import Partition_Command
from util.Ungap_Command import Ungap_Command
from viz.Visualize_OTU_Heatmap_Command import Visualize_OTU_Heatmap_Command
from viz.Visualize_OTU_Sample_Composition_Command import Visualize_OTU_Sample_Composition_Command

# Program Version
version = "0.01"
# Time format for printing
FORMAT = "%(asctime)s  %(message)s"
# Date format for printing
DATEFMT = "%m/%d %H:%M:%S"


# unused leters: h,k
# TODO parser args should probably live next to their program strings (I.E. in Program Runner or another class)
def main(argv):
    """Parses command line args, builds an argparse.ArgumentParser, and runs the chosen command.
        Otherwise, prints usage.

    :param argv: Command line arguments as a list of strings
    """
    parser = argparse.ArgumentParser(description="arms description", epilog="arms long description")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    parser.add_argument("--verbose", help="Increase output verbosity")
    parser.add_argument('-t', '--processes', type=int, default=1, help="The maximum number of processes to spawn.")
    parser.add_argument('--dry_run', default=False)
    parser.add_argument('--debugtest', default=False)
    parser.add_argument('-y', '--extraargstring', default="", required=False, help="Auxillary command line arguments \
                            to append to the end of the command.")

    subparsers = parser.add_subparsers(dest='action', help='Available commands')

    # ====================================
    # ==  Fix reads with Baye's hammer  ==
    # ====================================
    # "SPADES_PRECLEAN":  --only-error-correction -o %s -1 %s -2 %s"
    parser_preclean = subparsers.add_parser('preclean', description="Given a pair of left and right fasta/fastq reads, \
                            or a pair of folder containing the left and right fasta/fastq files, fixes short errors \
                            in the reads using Baye's Hammer.  Forwards reads filegroups should end in \
                            '_forward.<ext>' or '_R1.<ext>'.  Reverse reads filegroups should end in '_reverse.<ext>' \
                            or '_R2.<ext>' (where <ext> is the file extension).")
    parser_preclean.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads file or folder.")
    parser_preclean.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads file or folder.")
    parser_preclean.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_preclean.add_argument('-p', '--program', required=False, default="bayeshammer", help="Indicates which \
                            program to use.  Choices are: 'bayeshammer'.  Default: 'bayeshammer'.")
    parser_preclean.add_argument('-j', '--bayesthreads', type=int, required=False, default=1, help="The number of \
                            threads to use per process (default is 1")
    parser_preclean.set_defaults(command=Preclean_Command)

    # =================================
    # ==  Assemble Reads using pear  ==
    # =================================
    # "pear": programPaths["PEAR"] + " -f \"%s\" -r \"%s\" -o \"%s\" -j %d "
    parser_assemble = subparsers.add_parser('assemble', description="Given a pair of left and right fasta/fastq reads, \
                            or a pair of folder containing the left and right fasta/fastq files, assembles the left \
                            and right read files into contiguous sequence files.  Forwards reads filegroups should \
                            end in '_forward.<ext>' or '_R1.<ext>'.  Reverse reads filegroups should end in \
                            '_reverse.<ext>' or '_R2.<ext>'.  (where <ext> is the file extension).")
    parser_assemble.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads file or folder.")
    parser_assemble.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads file or folder.")
    parser_assemble.add_argument('-n', '--name', required=True, help="Assembled File Prefix.")
    parser_assemble.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_assemble.add_argument('-p', '--program', required=False, default="pear", help="Indicates which \
                            program to use.  Choices are: 'pear'.  Default: 'pear'.")
    parser_assemble.add_argument('-j', '--pearthreads', type=int, required=False, default=1, help="Pear: The number of \
                            threads to use per process (default is 1")
    parser_assemble.set_defaults(command=Assemble_Command)

    # ====================================
    # ==   Split by barcode with FastX  ==
    # ====================================
    # "barcode.splitter": "cat \"%s\" | " + programPaths["FASTX"] + "fastx_barcode_splitter.pl  --bcfile \"%s\" \
    #                                    -prefix \"%s\" --suffix .fastq --bol --mismatches 1",
    parser_demux = subparsers.add_parser('demux_samples', description="Given a single barcodes file, and a \
                            fasta/fastq file or a folder containing fasta/fastq files, splits each input fasta/fastq \
                            file into separate sample files (based on each sequence's barcode), for each barcode in \
                            the barcodes file.  Sample files of size zero are ignored.")
    parser_demux.add_argument('-i', '--input_f', required=True, help="Input fasta/fastq file or folder.")
    parser_demux.add_argument('-b', '--barcodes', required=True,
                              help="Tab delimted files of barcodes and corresponding samples.")
    parser_demux.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_demux.add_argument('-p', '--program', required=False, default="fastx", help="Indicates which \
                            program to use.  Choices are: 'fastx'.  Default: 'fastx'.")
    parser_demux.set_defaults(command=Demux_Command)

    # ==================================================
    # ==  Rename reads serially with renameSequences  ==
    # ==================================================
    # renameSequences(input, output)
    parser_rename = subparsers.add_parser('rename', description="Given a fasta/fastq file or directory of fasta/fastq \
                            files, serially regroups each sequence in each file serially, with the filename as a \
                            prefix.")
    parser_rename.add_argument('-i', '--input_f', required=True, help="Input fasta/fastq file or folder.")
    parser_rename.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_rename.add_argument('-f', '--filetype', required=True, help="The filetype of the input files.  Either \
                            'fasta' or 'fastq'.")
    parser_rename.add_argument('-c', '--clip', required=False, default=True, help="Set True if input file groups \
                            contain trailing demux_seqs identifiers.  e.g. True if file name contains '_0', '_1', \
                            '_2', etc..  Default: True.")
    parser_rename.add_argument('-p', '--program', required=False, default="chewbacca", help="Indicates which \
                            program to use.  Choices are: 'chewbacca'.  Default: 'chewbacca'.")
    parser_rename.set_defaults(command=Rename_Command)

    # =================================================
    # ==  Trims barcodes and adapters using flexbar  ==
    # =================================================
    # "flexbar":  "flexbar -r \"%s\" -t \"%s\" -ae \"%s\" -a \"%s\""
    parser_trim = subparsers.add_parser('trim_adapters', description="Given a single barcodes file, a single adapters \
                            file, and a single fasta/fastq file or a folder containing fasta/fastq files, removes the \
                            specified adapters (and preceeding barcodes) from all sequences in the given fasta/fastq \
                            files.")
    parser_trim.add_argument('-i', '--input_f', required=True, help="Input fasta/fastq file or folder.")
    parser_trim.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_trim.add_argument('-a', '--adapters', required=True, help="Forwards Adapters file.")
    parser_trim.add_argument('-arc', '--adaptersrc', required=True, help="Reverse Complimented Adapters file.")
    parser_trim.add_argument('-p', '--program', required=False, default="flexbar", help="Indicates which \
                            program to use.  Choices are: 'flexbar'.  Default: 'flexbar'.")
    parser_trim.add_argument('-u', '--allowedns', required=False, default=0, type=int, help=" Flexbar: The number of \
                            unknown 'N' bases a sequence is allowed before being thrown out.  Default: 0.")
    parser_trim.set_defaults(command=Clean_Adapters_Command)

    # ====================================
    # ==  Clean Reads with Low Quality  ==
    # ====================================
    # Clean low-quality reads with trimmomatic
    # "trimomatic":       "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
    # -phred33 input output_cleaned.fastq SLIDINGWINDOW:%windowsize:%minAvgQuality MINLEN:%minLen"
    parser_trimmomatic = subparsers.add_parser('clean_seqs', description="Given a single fastq file or a folder \
                            containing fastq files, trims away areas of low read quality (specified by -q) in each \
                            sequence.  Then, the longest remaining segment (of at least length -m) is recorded as the \
                            new sequence.")
    parser_trimmomatic.add_argument('-i', '--input_f', required=True, help="Input fastq file or folder")
    parser_trimmomatic.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_trimmomatic.add_argument('-p', '--program', required=False, default="trimmomatic", help="Indicates which \
                            program to use.  Choices are: 'trimmomatic'.  Default: 'trimmomatic'.")
    parser_trimmomatic.add_argument('-m', '--minlen', type=int, default=200,
                                    help="Trimmomatic: Minimum length for cleaned sequences")
    parser_trimmomatic.add_argument('-w', '--windowSize', type=int, default=5,
                                    help="Trimmomatic: Size of the sliding window")
    parser_trimmomatic.add_argument('-q', '--quality', type=int, default=25,
                                    help="Trimmomatic: Minimum average quality for items in the sliding window")
    parser_trimmomatic.set_defaults(command=Clean_Quality_Command)

    # ==========================================
    # ==  Dereplicate sequences with vsearch  ==
    # ==========================================
    # " vsearch --threads %d --derep_fulllength %s --sizeout --fasta_width 0 --output %s -uc %s",
    parser_derep = subparsers.add_parser('dereplicate_fasta', description="Given an fasta file or folder, removes \
                            identical duplicates and subsequences WITHIN EACH FILE, keeping the longest sequence, and \
                            regroups it with the number of duplicates as '<longest_sequence_name>_<duplicate count>'.")
    parser_derep.add_argument('-i', '--input_f', required=True, help="Input fasta file or folder of fasta files.")
    parser_derep.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_derep.add_argument('-j', '--threads', required=False, type=int, default=2, help="Number of threads to use \
                            per query process.")
    parser_derep.add_argument('-g', '--groupsfile', required=False, help="A .groups file to update.  If no .groups \
                            file is provided, then sequences are assumed to be singletons.")
    parser_derep.add_argument('-p', '--program', required=False, default="vsearch", help="Indicates which \
                            program to use.  Choices are: 'vsearch'.  Default: 'vsearch'.")
    parser_derep.add_argument('-s', '--stripcounts', required=False, type=bool, default=False, help="If included, \
                            strip counts from sequence groups before clustering.  This allows for the recognition of \
                            sequence groups that are annotated with dereplication counts.")
    parser_derep.set_defaults(command=Dereplicate_Command)

    # ============================================
    # ==  Partition fastas with splitKperFasta  ==
    # ============================================
    # splitK(inputFasta, prefix, nbSeqsPerFile, filetype):
    parser_split = subparsers.add_parser('partition', description="Given a fasta/fastq file, or folder containing \
                            fasta/fastq files, splits each file into a set of sub files labeled \
                            <filename>_part_x.<ext>, where each subfile contains at most <--chunksize> sequences. \
                            (Where <ext> is the file extension of the original file.")
    parser_split.add_argument('-i', '--input_f', required=True, help="Input fasta/fastq file or folder.")
    parser_split.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_split.add_argument('-c', '--chunksize', type=int, required=True, help="Chunksize. (Max number of sequences \
                            per file).")
    parser_split.add_argument('-f', '--filetype', required=True, help="The filetype of the input files.  Either \
                            'fasta' or 'fastq'.")
    parser_split.add_argument('-p', '--program', required=False, default="chewbacca", help="Indicates which \
                            program to use.  Choices are: 'chewbacca'.  Default: 'chewbacca'.")
    parser_split.set_defaults(command=Partition_Command)

    # =======================
    # ==  Merge the files  ==
    # =======================
    # merge_files(input_file_list, output_file)
    parser_cat = subparsers.add_parser('merge_files', description="Given a folder containing any type of file, \
                            concatenates all files in that folder together, regardless of their filetype.  Note: Any \
                            file headers/footers will be written as well.")
    parser_cat.add_argument('-i', '--input_f', required=True, help="Input directory containing files to concatenate.")
    parser_cat.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_cat.add_argument('-n', '--name', required=True, help="Prefix name for the merged file.")
    parser_cat.add_argument('-f', '--fileext', required=True, help="File extension for the output file.")
    parser_cat.add_argument('-p', '--program', required=False, default="chewbacca", help="Indicates which \
                            program to use.  Choices are: 'chewbacca'.  Default: 'chewbacca'.")
    parser_cat.set_defaults(command=Merge_Command)

    # =======================
    # ==  ungap the files  ==
    # =======================
    # ungap(file_to_clean, output_file_name, gap_char, file_type)
    parser_ungap = subparsers.add_parser('ungap_fasta', description="Given a fasta/fastq file or a folder containing \
                            fasta/fastq files, removes alignment characters (specified by -g).")
    parser_ungap.add_argument('-i', '--input_f', required=True, help="Input directory containing files to concatenate.")
    parser_ungap.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_ungap.add_argument('-f', '--fileext', required=True, help="File extension for the output file.  Either \
                            'fasta', or 'fastq'.")
    parser_ungap.add_argument('-g', '--gapchars', required=True, help="A string of one or more characters to remove \
                            from the sequences (but not sequence groups) in the input files.")
    parser_ungap.add_argument('-p', '--program', required=False, default="chewbacca", help="Indicates which \
                            program to use.  Choices are: 'chewbacca'.  Default: 'chewbacca'.")
    parser_ungap.set_defaults(command=Ungap_Command)

    # ==============================
    # ==  Align Reads with MACSE  ==
    # ==============================
    # "macse_align":      "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
    #                                \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
    #                                -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
    #                                -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",
    parser_align = subparsers.add_parser('macseAlign', description="Aligns sequences using the 'enrichAlignment' \
                            program in Macse.  Removes chimeras, invalid sequences, and orients sequences in the \
                            forward read direction.")
    parser_align.add_argument('-i', '--input_f', required=True, help="Input fasta")
    parser_align.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_align.add_argument('-d', '--db', required=True, help="Database against which to align and filter reads")
    parser_align.add_argument('-p', '--program', required=False, default="macse", help="Indicates which \
                            program to use.  Choices are: 'macse'.  Default: 'macse'.")
    parser_align.set_defaults(command=Clean_Deep_Command)

    # =====================================
    # == Clean Aligned Reads with MACSE  ==
    # =====================================
    # "macse_format":     "java -jar " + programPaths["MACSE"] + "  -prog exportAlignment -align \"%s\" \
    #                                -charForRemainingFS - -gc_def 5 -out_AA \"%s\" -out_NT \"%s\" -statFile \"%s\"",
    parser_macse_clean = subparsers.add_parser('macseClean', description="Removes gap characters inserted by the macse \
                            aligner.")
    parser_macse_clean.add_argument('-i', '--input_f', required=True, help="Input fasta")
    parser_macse_clean.add_argument('-s', '--samplesdir', required=True, help="Samples dir")
    parser_macse_clean.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_macse_clean.add_argument('-d', '--db', required=True, help="Database against which to align and filter \
                            reads")
    parser_macse_clean.add_argument('-p', '--program', required=False, default="macse", help="Indicates which \
                             program to use.  Choices are: 'macse'.  Default: 'macse'.")
    parser_macse_clean.set_defaults(command=Clean_Deep_Repair_Command)

    # ============================================
    # ==  Cluster using CROP, SWARM, OR VSEARCH ==
    # ============================================
    parser_cluster = subparsers.add_parser('cluster_seqs', description="Given a fasta file or a folder containing \
                            fasta files, performs clustering on each input file individually.  Outputs a .groups file \
                            listing unique representative sequences from each cluster, and the groups of the\
                            sequences they represent.  Also outputs a <input_file_name>_seeds.fasta file of the unique \
                            representatives.")
    parser_cluster.add_argument('-i', '--input_f', required=True, help="Input fasta file/folder")
    parser_cluster.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_cluster.add_argument('-p', '--program', required=False, default="swarm", help="One of 'crop', 'swarm', or \
                            'vsearch' to indicate which clustering program to use.  Default: 'swarm'.")
    parser_cluster.add_argument('-g', '--groupsfile', required=False, help="A .groups file to update.")
    parser_cluster.add_argument('-s', '--stripcounts', required=False, type=bool, default=True, help="If True, strip \
                            counts from sequence groups before clustering.  This allows for the recognition of \
                            group names in cases where dereplication counts have been added to group names.")
    # CROP options
    parser_cluster.add_argument('-z', '--blocksize', required=False, type=int, default=500, help="CROP only: Size of \
                            blocks to be used for all rounds (if -b is specified, then -z will not affect the first \
                            round.  For data set with different average sequence length, this parameter should be \
                            tuned such that it won't take too long for each block to do pariwise alignment.  Hint for \
                            choosing z: z*L<150,000, where L is the average length of the sequences.  Default: 500.")
    parser_cluster.add_argument('-b', '--blockcount', required=False, type=int, help="CROP only: The size of blocks in \
                            the first round of clustering. Hint of choosing -b: Each block in the first round should \
                            contain about 50 sequences.  i.e. b=N/50, where N is the number of input sequences.  \
                            Default: # input sequences / z.")
    parser_cluster.add_argument('-e', '--maxmcmc', required=False, type=int, default=2000, help="CROP only: This \
                            parameter specifies the number of iterations of MCMC. Default value is 2000. Increase this \
                            value to enhance accuracy (recommended value is at least 10*block size).")
    parser_cluster.add_argument('-c', '--clustpct', required=False, default="g", help="CROP only: The minimum \
                            similarity threshold for clustering.  Either 'g' for 95% or 's' for 97%.  Default: 'g'.")
    parser_cluster.add_argument('-m', '--maxsm', required=False, type=int, default=20, help="CROP only: This parameter \
                            specifies the maximum number of 'split and merge' process to run. Default value is 20, \
                            which is also the maximum allowed.")
    parser_cluster.add_argument('-r', '--rare', required=False, type=int, default=2, help="CROP only: The maximum \
                            cluster size allowed to be classified as 'rare'. Clusters are defined as either \
                            'abundant' or 'rare'. 'Abundant' clusters will be clustered first, then the 'rare' \
                            clusters are mapped to the 'abundant' clusters.  Finally, 'rare' clusters which cannot be \
                            mapped will be clustered separately. e.g. If r=5, the clusters with size <=5 will be \
                            considered 'rare' in above procedure. and r=0 will yield the best accuracy. If you \
                            believe your data is not too diverse to be handled, then r=0 will be the best choice. \
                            Default: 2.")
    # Vsearch options
    parser_cluster.add_argument('-v', '--idpct', required=False, type=float, default=.95, help="VSEARCH only: % match \
                            required for clustering.  Real number in the range (0,1]. Default: 0.95")
    parser_cluster.set_defaults(command=Cluster_Command)

    # ================================================
    # ==  Closed Ref Picking with fasta references  ==
    # ================================================
    # --usearch_global  ../9_p_uchime/seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
    # --userfields query+target+id+alnlen+qcov --userout out  --alnout alnout.txt
    parser_query_fasta = subparsers.add_parser('query_fasta', description="Given a fasta file, or folder containing \
                            fasta files, aligns and identifies OTUs against a curated fasta file.")
    parser_query_fasta.add_argument('-i', '--input_f', required=True, help="Input file/folder with fasta files")
    parser_query_fasta.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_query_fasta.add_argument('-r', '--referencefasta', required=True, help="Filepath to the curated fasta file \
                            to use as a reference.")
    parser_query_fasta.add_argument('-x', '--taxinfo', required=True, help="Filepath to a two-column, tab-delimited \
                            file mapping a sequence's fasta id (in the referencefasta file) to a taxonomic \
                            identification.")
    parser_query_fasta.add_argument('-s', '--simmilarity', required=False, default=97, type=int, help="Minimum %  \
                            simmilarity (integer between 0 and 100) between query and reference sequences required for \
                            positive identification. Default: 97")
    parser_query_fasta.add_argument('-p', '--program', required=False, default="vsearch", help="Indicates which \
                            program to use.  Choices are: 'vsearch'.  Default: 'vsearch'.")
    parser_query_fasta.add_argument('-c', '--coverage', required=False, default=85, type=int, help="Minimum % coverage \
                            (integer between 0 and 100) required query and reference sequences required for positive \
                            identification. Default: 85")
    parser_query_fasta.add_argument('-j', '--threads', required=False, type=int, default=2,
                                    help="Number of threads to use per query process.")
    parser_query_fasta.set_defaults(command=Query_OTU_Fasta_Command)

    # ====================================
    # ==  Closed Ref Picking with NCBI  ==
    # ====================================
    # --usearch_global ../9_p_uchime/seeds.pick.fasta  --db /home/mahdi/refs/COI_DOWNLOADED/COI.fasta -id 0.9 \
    #          --userfields query+target+id+alnlen+qcov --userout out --alnout alnout.txt
    # --userfields query+target+id+alnlen+qcov
    parser_query_db = subparsers.add_parser('query_db', description="Given a fasta file, or folder containing \
                            fasta files, aligns and identifies OTUs against a curated database.")
    parser_query_db.add_argument('-i', '--input_f', required=True, help="Input Fasta File to clean")
    parser_query_db.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved")
    parser_query_db.add_argument('-r', '--referencefasta', required=True, help="Filepath to the curated fasta file \
                            to use as a reference.")
    parser_query_db.add_argument('-d', '--db', required=True, help="Filepath to the curated fasta file \
                            to use as a reference.")
    parser_query_db.add_argument('-s', '--simmilarity', required=False, default=97, type=int, help="Minimum %  \
                                simmilarity (integer between 0 and 100) between query and reference sequences \
                                required for positive identification. Default: 97")
    parser_query_db.add_argument('-c', '--coverage', required=False, default=85, type=int, help="Minimum % coverage \
                                (integer between 0 and 100) required query and reference sequences required for \
                                positive identification. Default: 85")
    parser_query_db.add_argument('-j', '--threads', required=False, type=int, default=2,
                                 help="Number of threads to use per query process.")
    parser_query_db.add_argument('-p', '--program', required=False, default="vsearch", help="Indicates which \
                            program to use.  Choices are: 'vsearch'.  Default: 'vsearch'.")
    parser_query_db.set_defaults(command=Query_OTU_DB_Command)

    # =======================
    # ==  Build OTU table  ==
    # =======================
    parser_build_matrix = subparsers.add_parser('build_matrix', description="Given a single barcodes file with all possible \
                            sample groups, a list of the latest groups file(s), and a list of initial samples files \
                            (mapping each original, undereplicated sequence to its sample name), builds an OTU \
                            table..  ")
    parser_build_matrix.add_argument('-s', '--samples', required=True, help="Input samples file or folder.")
    parser_build_matrix.add_argument('-g', '--groups', required=True, help="Input groups file or folder.")
    parser_build_matrix.add_argument('-b', '--barcodes', required=True, help="Input barcodes file.")
    parser_build_matrix.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_build_matrix.add_argument('-p', '--program', required=False, default="chewbacca", help="Indicates which \
                            program to use.  Choices are: 'chewbacca'.  Default: 'chewbacca'.")
    parser_build_matrix.set_defaults(command=Build_OTU_Table_Command)

    # =======================
    # ==  Annotate Matrix  ==
    # =======================
    parser_annotate_matrix = subparsers.add_parser('annotate_matrix', description="Given a tabular file mapping \
                            sequence IDs to taxonomic groups, and an OTU matrix, regroups the identifiable \
                            sequence IDs with taxonomic groups.")
    parser_annotate_matrix.add_argument('-i', '--input_f', required=True,
                                        help="Input matrix file or folder of matrix files.")
    parser_annotate_matrix.add_argument('-a', '--annotation', required=True,
                                        help="File mapping sequence IDs to taxonomic groups.")
    parser_annotate_matrix.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_annotate_matrix.add_argument('-p', '--program', required=False, default="chewbacca", help="Indicates which \
                            program to use.  Choices are: 'chewbacca'.  Default: 'chewbacca'.")
    parser_annotate_matrix.set_defaults(command=Annotate_OTU_Table_Command)

    # ==============================
    # ==  Convert fastq to fasta  ==
    # ==============================
    # translateFastqToFasta(inputFastQ, outputFasta):
    parser_to_fasta = subparsers.add_parser('convert_fastq_to_fasta', description="Converts a fastq file to a fasta \
                            file.")
    parser_to_fasta.add_argument('-i', '--input_f', required=True, help="Fastq file or folder containing fastq files to \
                            translate")
    parser_to_fasta.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_to_fasta.add_argument('-p', '--program', required=False, default="chewbacca", help="Indicates which \
                            program to use.  Choices are: 'chewbacca'.  Default: 'chewbacca'.")
    parser_to_fasta.set_defaults(command=Convert_Fastq_Fasta_Command)

    # ===========================================
    # == Visualize Samples as OTU Composition  ==
    # ===========================================
    parser_viz_otu_comp = subparsers.add_parser('visualize_otu_sample_composition', description="Creates a stacked \
                            barchart showing the OTU composition in each sample.")
    parser_viz_otu_comp.add_argument('-i', '--input_f', required=True, help="Filepath to the OTU table to visualize.  \
                            Input OTU tables should start with a tab-delimited header row of samplenames, prefixed by \
                            the word \"OTU\".  Each subsequent line should be a tab-delimited line listing an OTU, \
                            followed by its abundance in each Sample.\ne.g.:\n\
                            OTU <Sample_name1> <Sample_name2> <Sample_name3>...\n\
                            otu_name1 0 1 5...\n\
                            otu_name2 1 2 0...\n\
                            otu_name3 3 1 1...\n")
    parser_viz_otu_comp.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_viz_otu_comp.add_argument('-p', '--program', required=False, default="chewbacca", help="Indicates which \
                            program to use.  Choices are: 'chewbacca'.  Default: 'chewbacca'.")
    group = parser_viz_otu_comp.add_mutually_exclusive_group(required=False)
    group.add_argument('-m', '--pct', help="Real number x in the range (0,1] indicating that the top x% of OTU names \
                           (sorted by abundance) should be included in the graph.", type=float)
    group.add_argument('-n', '--names', help="Filepath to a list of OTU names to include.  File should be formatted as \
                            a series of lines where each line contains just an OTU name.")
    group.add_argument('-k', '--count', help="Positive integer x indicating that x OTU names (sorted highest \
                            abundance) should be included in the graph.", type=int)
    group.set_defaults(command=Visualize_OTU_Sample_Composition_Command)

    # ===========================
    # == Visualize OTU Heatmap ==
    # ===========================
    parser_viz_otu_heatmap = subparsers.add_parser('visualize_otu_heatmap', description="Creates a heatmap showing the \
                            abundance of each OTU at each sample site.")
    parser_viz_otu_heatmap.add_argument('-i', '--input_f', required=True, help="Filepath to the OTU table to visualize.  \
                            Input OTU tables should start with a tab-delimited header row of samplenames, prefixed by \
                            the word \"OTU\".  Each subsequent line should be a tab-delimited line listing an OTU, \
                            followed by its abundance in each Sample.\ne.g.:\n\
                            OTU <Sample_name1> <Sample_name2> <Sample_name3>...\n\
                            otu_name1 0 1 5...\n\
                            otu_name2 1 2 0...\n\
                            otu_name3 3 1 1...\n")
    parser_viz_otu_heatmap.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    parser_viz_otu_heatmap.add_argument('-p', '--program', required=False, default="chewbacca", help="Indicates which \
                            program to use.  Choices are: 'chewbacca'.  Default: 'chewbacca'.")
    group = parser_viz_otu_heatmap.add_mutually_exclusive_group(required=False)
    group.add_argument('-m', '--pct', help="Real number x in the range (0,1] indicating that the top x% of OTU names \
                           (sorted by abundance) should be included in the graph.", type=float)
    group.add_argument('-n', '--names', help="Filepath to a list of OTU names to include.  File should be formatted as \
                            a series of lines where each line contains just an OTU name.")
    group.add_argument('-k', '--count', help="Positive integer x indicating that x OTU names (sorted highest \
                            abundance) should be included in the graph.", type=int)
    parser_viz_otu_heatmap.set_defaults(command=Visualize_OTU_Heatmap_Command)

    # =======================================
    # == Parse args and call default func  ==
    # =======================================
    args, unknown = parser.parse_known_args()
    if unknown:
        print "\nIgnoring unknown args: " + ', '.join(['%s'] * len(unknown)) % tuple(unknown)
    if args.verbose:
        logging.basicConfig(format=FORMAT, level=logging.DEBUG, datefmt=DATEFMT)
    else:
        logging.basicConfig(format=FORMAT, level=logging.ERROR, datefmt=DATEFMT)

    printVerbose.VERBOSE = (args.verbose is not None)
    logging.debug("Initial ARGS are: %s", args)
    print("\t\t")
    dryRun = args.dry_run
    signal.signal(signal.SIGTSTP, signal.SIG_IGN)
    makeDirOrdie(args.outdir)
    args.command(args).execute_command()


if __name__ == "__main__":
    main(sys.argv)
