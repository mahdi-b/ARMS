from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *
from util.ungap import remove_gap_chars

from classes.Helpers import *


class Ungap_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"


    def execute_program(self):
        args = self.args
        self.ungap_chewbacca(args.input_f, args.outdir, args.gapchars, args.fileext, args.processes)


    def ungap_chewbacca(self, input_f, outdir, gapchars, file_ext, processes):
        """Ungaps a character using Bio python.

            :param input_f: Filepath to input file or folder to ungap.
            :param outdir: Filepath to the output directory where ungapped files should be written.
            :param gapchars: A string containing the gap characters to remove.
            :param file_ext: Either 'fasta' or 'fastq'.
            :param processes: The number of threads to use to ungap the input fileset.
        """

        makeDirOrdie(outdir)
        inputs = getInputFiles(input_f, "*.fasta")
        debugPrintInputInfo(inputs, "ungap.")
        pool = init_pool(min(len(inputs), processes))
        printVerbose("Removing all '%s' from sequences..." % gapchars)
        # ungap(file_to_clean, output_file_name, gap_char, file_type):
        parallel(runPythonInstance, [(remove_gap_chars, input_, "%s/%s_cleaned.%s" % (outdir, strip_ixes(input_), 'fasta'),
                                      gapchars, file_ext) for input_ in inputs], pool)
        printVerbose("Done removing.")
        cleanup_pool(pool)
