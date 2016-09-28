from classes.Helpers import *
from classes.ProgramRunner import *
from clean_gap_chars_chewbacca import remove_gap_chars

def ungap_main(input_f, outdir, gapchars, file_ext, program, threads, aux_params):
    """Removes a gap character from sequences in a fasta file.  Useful for removing characters from an alignment file.

    :param input_f: Filepath to input file or folder to ungap.
    :param outdir: Filepath to the output directory where ungapped files should be written.
    :param gapchars: A string containing the gap characters to remove.
    :param file_ext: Either 'fasta' or 'fastq'.
    :param program: The name of the program to use.  Choices are ["chewbacca"].
    :param threads: The number of threads to use to ungap the input fileset.
    :param aux_params: A dictionary of program-specific named-parameters.
    """
    makeDirOrdie(outdir)
    inputs = getInputFiles(input_f, "*.fasta")
    debugPrintInputInfo(inputs, "ungap.")
    pool = init_pool(min(len(inputs), threads))

    if program == "chewbacca":
        ungap_chewbacca(inputs, outdir, gapchars, threads, file_ext, pool)
    cleanup_pool(pool)


def ungap_chewbacca(inputs, outdir, gapchars, file_ext, pool):
    """Ungaps a character using Bio python.

    :param inputs: A list of filepaths to input files to ungap.
    :param outdir: Filepath to the output directory where ungapped files should be written.
    :param gapchars: A string containing the gap characters to remove.
    :param file_ext: Either 'fasta' or 'fastq'.
    :param pool: A fully-initalized multiprocessing.Pool object.
    """
    printVerbose("Removing all '%s' from sequences..." % gapchars)
    # ungap(file_to_clean, output_file_name, gap_char, file_type):
    parallel(runPythonInstance, [(remove_gap_chars, input_, "%s/%s_cleaned.%s" % (outdir, strip_ixes(input_), 'fasta'),
                                  gapchars, file_ext) for input_ in inputs], pool)
    printVerbose("Done removing.")
