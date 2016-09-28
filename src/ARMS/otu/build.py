from classes.Helpers import *
from buildOTUtable import buildOTUtable


def build_otu_main(outdir, groups_file, samples_file, program, barcodes_file, aux_params):
    """Builds the unannotated OTU table using custom chewbacca script.

    :param outdir: The directory where the matrix should be written.
    :param groups_file: A .groups file containing the OTU names and their consituent/replicant sequences.
    :param samples_file: A .samples file containing the samples that each sequence in the .groups file belongs to.
    :param barcodes_file: A .barcodes file listing all sample names.
    :param aux_params: A dictionary of program-specifc named-parameters.
    """
    makeDirOrdie(outdir)
    if program == "chewbacca":
        build_otu_chewbacca(outdir, groups_file, samples_file, barcodes_file)


def build_otu_chewbacca(outdir, groups_file, samples_file, barcodes_file):
    """Builds the unannotated OTU table using custom chewbacca script.

    :param outdir: The directory where the matrix should be written.
    :param groups_file: A .groups file containing the OTU names and their consituent/replicant sequences.
    :param samples_file: A .samples file containing the samples that each sequence in the .groups file belongs to.
    :param barcodes_file: A .barcodes file listing all sample names.
    """
    groups = getInputFiles(groups_file)
    debugPrintInputInfo(groups, "read.")
    samples = getInputFiles(samples_file)
    debugPrintInputInfo(samples, "read.")
    barcodes = getInputFiles(barcodes_file)
    debugPrintInputInfo(barcodes, "read.")
    printVerbose("Building matrix...")
    buildOTUtable(groups, samples, barcodes[0], "%s/%s.txt" % (outdir, "matrix"))
    printVerbose("Done building.")