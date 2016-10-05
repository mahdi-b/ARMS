from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *

from classes.Helpers import *


class Demux_Program_Fastx(ChewbaccaProgram):
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
    """
    name = "fastx"

    def execute_program(self):
        args = self.args
        self.demux_fastx(args.input_f, args.barcodes, args.outdir, args.processes, args.extraargstring)

    def demux_fastx(self, input_f, barcodes, outdir, processes, extraargstring):
        """Demuxes using FAST X BARCODE SPLITTER.

        :param input_f: File path to input file or folder of input files.
        :param barcodes: File path to input barcodes file.
        :param outdir: Filepath to output directory.
        :param processes: Number of processes to use to demux input fileset.
        :param extraargstring: Advanced program parameter string.
        """
        makeDirOrdie(outdir)
        # Get input files
        files_to_split = getInputFiles(input_f)
        # Assign the files shard numbers
        file_id = range(len(files_to_split))
        file_id_pairs = zip(files_to_split, file_id)
        debugPrintInputInfo(files_to_split, "demux")
        pool = init_pool(min(len(file_id_pairs), processes))

        printVerbose("Demuxing sequences...")
        parallel(runProgramRunnerInstance,
                 [ProgramRunner(ProgramRunnerCommands.DEMUX_FASTX,
                                [input_, barcodes, "%s/" % outdir, "_%d_splitOut.fastq" % id_],
                                {"exists": [input_, barcodes]}, extraargstring)
                  for input_, id_ in file_id_pairs], pool)
        printVerbose("Demuxed sequences.")

        # gather output files and move them to their final destination
        output_files = getInputFiles("*_splitOut_", ignore_empty_files=False)
        bulk_move_to_dir(output_files, outdir)

        # Grab all the auxillary files
        aux_files = getInputFiles(outdir, "unmatched_*", ignore_empty_files=False)
        # make aux dir for extraneous files and move them there
        bulk_move_to_dir(aux_files, makeAuxDir(outdir))
        cleanup_pool(pool)
