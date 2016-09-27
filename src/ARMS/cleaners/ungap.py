from classes.Helpers import *
from classes.ProgramRunner import *
from remove_gaps import remove_gap_chars

def ungap_fasta(input_f, outdir, gapchars, file_ext, program, threads, aux_params):
    """Removes a gap character from sequences in a fasta file.  Useful for removing characters from an alignment file.
       :param args: An argparse object with the following parameters:
                       input        Fasta files to ungap.
                       outdir       Directory where outputs will be saved.
                       gapchar      A gap character to remove from the sequences.
                       fileext      Either 'fastq' or 'fasta'.
       """
    makeDirOrdie(outdir)
    inputs = getInputFiles(input_f, "*.fasta")
    debugPrintInputInfo(inputs, "ungap.")
    pool = init_pool(min(len(inputs), threads))

    if program == "chewbacca":
        ungap_chewbacca(inputs, outdir, gapchars, threads, file_ext, pool)
    cleanup_pool(pool)


def ungap_chewbacca(inputs, outdir, gapchars, threds, file_ext, pool):
    printVerbose("Removing all '%s' from sequences..." % gapchars)
    # ungap(file_to_clean, output_file_name, gap_char, file_type):
    parallel(runPythonInstance, [(remove_gap_chars, input_, "%s/%s_cleaned.%s" % (outdir, strip_ixes(input_), 'fasta'),
                                  gapchars, file_ext) for input_ in inputs], pool)
    printVerbose("Done removing.")
