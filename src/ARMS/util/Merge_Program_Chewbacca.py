from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *
from util.merge import merge_files

from classes.Helpers import *


class Merge_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"


    def execute_program(self):
        args = self.args
        self.merge_chewbacca(args.input_f, args.outdir, args.name, args.fileext)



    def merge_chewbacca(self, input_f, outdir, output_filename, output_fileext):
        """Merges files together into a new output file.

        :param input_f: Filepath to a directory of input files.
        :param outdir: Filepath to the output folder.
        :param program: The program to use to merge files.  Choices are ["chewbacca"]. Default: "chewbacca".
        :param output_filename: The filename of the output file, without an extension.
        :param output_fileext: The file extension of the output file.
        :param aux_params: A dictionary of program-specific named-parameters.
        """
        makeDirOrdie(outdir)
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "merged")
        printVerbose("Merging files.")
        output_file = "%s/%s_MERGED.%s" % (outdir, output_filename, output_fileext)
        merge_files(inputs, output_file)
        printVerbose("Done merging.")
        logging.debug("Merged %d files into %s" % (len(inputs), output_file))