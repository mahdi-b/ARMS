from classes.Helpers import *
from util.joinFiles import joinFiles

def merge_main(input_f, outdir, program, output_filename, output_fileext, aux_params):
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
    debugPrintInputInfo(inputs, "merged together.")
    if program == "chewbacca":
        joinFiles(inputs, "%s/%s_MERGED.%s" % (outdir, output_filename, output_fileext))