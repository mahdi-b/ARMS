from Bio import SeqIO
from classes.BufferedWriter import BufferedSeqWriter
from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import makeAuxDir, makeDirOrdie, getInputFiles, debugPrintInputInfo, init_pool, printVerbose,\
                             bulk_move_to_dir, cleanup_pool, run_parallel
from classes.PythonRunner import PythonRunner
from parse.parseBarcodesToDict import parse_barcodes_to_dict


class Demux_Program_Chewbacca(ChewbaccaProgram):
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
    name = "chewbacca"

    def execute_program(self):
        args = self.args
        self.demux_by_name(args.input_f, args.barcodes, args.outdir, "fastq", args.processes, args.extraargstring)

    def demux_by_name(self, input_f, barcodes, outdir, filetype, processes, extraargstring):
        """Demuxes using SeqIO.

        :param input_f: File path to input file or folder of input files.
        :param barcodes: File path to input barcodes file.
        :param outdir: Filepath to output directory.
        :param filetype: Either 'fasta' or 'fastq'.
        :param processes: Number of processes to use to demux input fileset.
        :param extraargstring: Advanced program parameter string.
        """
        makeDirOrdie(outdir)
        aux_dir = makeAuxDir(outdir)
        # Get input files
        files_to_split = getInputFiles(input_f)
        # Assign the files shard numbers
        file_id = range(len(files_to_split))
        file_id_pairs = zip(files_to_split, file_id)
        debugPrintInputInfo(files_to_split, "demux")
        pool = init_pool(min(len(file_id_pairs), processes))

        printVerbose("Demuxing sequences...")
        run_parallel([PythonRunner(split_on_name,
                                   [input_, barcodes, outdir, id_, filetype], {"exists": [input_]})
                        for input_, id_ in file_id_pairs], pool)


        # Grab all the auxillary files
        aux_files = getInputFiles(outdir, "unmatched_*", ignore_empty_files=False)
        # make aux dir for extraneous files and move them there
        bulk_move_to_dir(aux_files, aux_dir)
        cleanup_pool(pool)


def split_on_name(input_f, barcodes_file, outdir, id_, filetype):
    sample_names = parse_barcodes_to_dict(barcodes_file).keys()
    out_streams = {}
    for sample_name in sample_names:
        outfile =  "%s/%s_%s_splitOut.fastq" % ( outdir, sample_name, id_)
        out_streams[sample_name] = BufferedSeqWriter(outfile, filetype)

    SeqIO.parse(input_f, filetype)
    seq_dict = SeqIO.index(input_f, filetype)

    for name in seq_dict.keys():
        for sample_name in sample_names:
            if sample_name in name:
                out_streams[sample_name].write(seq_dict[name])

    for writer in out_streams.keys():
        out_streams[writer].flush()