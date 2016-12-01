import os
from classes.BufferedWriter import BufferedSeqWriter, BufferedFileWriter
from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.PythonRunner import PythonRunner
from classes.Helpers import getInputFiles, debugPrintInputInfo, init_pool, run_parallel, printVerbose, strip_ixes, \
    cleanup_pool, bulk_move_to_dir, makeAuxDir, makeDirOrdie, clip_count
from Bio import SeqIO


class Rename_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"

    def execute_program(self):
        args = self.args
        self.rename_chewbacca(args.input_f, args.outdir, args.filetype, args.clip, args.processes)

    def rename_chewbacca(self, input_f, outdir, filetype, clip, processes):
        """Renames sequences in a fasta/fastq file as <filename>_ID0, <filename>_ID1, <filename>_ID2, etc., where
            <filename> is the name of the fasta/fastq file without any extensions or chewbacca suffixes.

        :param input_f: Filepath to an input file or folder to rename.
        :param outdir: Filepath to the output directory.
        :param filetype: Either 'fasta' or 'fastq'.
        :param clip: If True, remove dereplication counts from sequence names before renaming.
        :param processes: The maximum number of processes to use.
        """

        # Gather input files
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "rename")
        pool = init_pool(min(len(inputs), processes))
        printVerbose("Renaming sequences...")
        # Run serialRename in run_parallel
        run_parallel([PythonRunner(serialRename,
                                   [input_,
                                    "%s/%s_renamed%s" % (outdir, strip_ixes(input_), os.path.splitext(input_)[1]),
                                    filetype, clip], {"exists": [input_]})
                      for input_ in inputs], pool)
        printVerbose("Done renaming sequences...")

        samples_dir = makeDirOrdie("%s_samples" % outdir)
        samples_files = getInputFiles(outdir, "*.samples", ignore_empty_files=False)
        bulk_move_to_dir(samples_files, samples_dir)

        aux_dir = makeAuxDir(outdir)
        aux_files = getInputFiles(outdir, "*.mapping", ignore_empty_files=False)
        bulk_move_to_dir(aux_files, aux_dir)

        cleanup_pool(pool)


def serialRename(input_file, output_fasta_filepath, file_type, clip=True):
    """Takes in a fasta file and outputs a new fasta with the sequences renamed.  Renaming convention is x.y.z<n> for
        x.y.z.fasta, where n is an integer in the range [0:n] where n is the position of the sequence in the input_file.
        Also writes a groups file, linking each sequence to its parent sample.
        e.g. The sequences in SiteX_SampleA.fasta are renamed:
                SiteX_SampleA_0, SiteX_SampleA_1, SiteX_SampleA_2, etc.
    :param input_file: Input fasta or fastq file.
    :param output_fasta_filepath:     Filepath for the output .samples file.
    :param file_type: "fasta" or "fastq"
    :param clip: True if filenames contain file_ID#s.  Will clip the IDs before renaming to get proper sequence names.
    """

    samples_file = "%s/%s_renamed.samples" % (os.path.dirname(output_fasta_filepath), strip_ixes(input_file))
    name_map_file = "%s/%s_renamed.mapping" % (os.path.dirname(output_fasta_filepath), strip_ixes(input_file))
    seq_prefix = strip_ixes(input_file)
    i = 0
    renamed_fasta = BufferedSeqWriter(output_fasta_filepath, file_type)
    mapping_file_output = BufferedFileWriter(name_map_file)
    samples_file_output = BufferedFileWriter(samples_file)

    for s in SeqIO.parse(input_file, file_type):
        i += 1

        # Store the old_name new_name mapping
        old_id = s.id
        s.id = "%s_ID%s" % (seq_prefix, i)
        mapping_file_output.write("%s\t%s" % (old_id, s.id))

        # Store the sequence-sample map
        if clip:
            sample_name = clip_count(seq_prefix)
        else:
            sample_name = seq_prefix
        samples_file_output.write("%s\t%s" % (s.id, sample_name))

        # Store the renamed sequence
        s.description = ""
        renamed_fasta.write(s)

    renamed_fasta.flush()
    mapping_file_output.flush()
    samples_file_output.flush()
