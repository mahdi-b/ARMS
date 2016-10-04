from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *

from classes.Helpers import *
from renameSerially import serialRename


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
        # Make the output directory, or abort if it already exists
        makeDirOrdie(outdir)

        # Gather input files
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "rename")
        pool = init_pool(min(len(inputs), processes))
        printVerbose("Renaming sequences...")
        # Run serialRename in parallel
        parallel(runPythonInstance,
                 [(serialRename,
                   input_, "%s/%s_renamed%s" % (outdir, strip_ixes(input_), os.path.splitext(input_)[1]),
                   filetype, clip) for input_ in inputs], pool)
        printVerbose("Done renaming sequences...")

        samples_dir = makeDirOrdie("%s_samples" % outdir)
        samples_files = getInputFiles(outdir, "*.samples", ignore_empty_files=False)
        bulk_move_to_dir(samples_files, samples_dir)
        cleanup_pool(pool)
