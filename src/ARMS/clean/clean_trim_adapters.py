from classes.Helpers import *
from clean_trim_adapters_flexbar import clean_trim_adapters_flexbar

def trim_adapters_main(input_f, adapters, adaptersrc, outdir, program, threads, aux_params):
    """Trim adapters (and preceeding barcodes) from sequences in input file(s).  Sequences should be in
        the following format: <BARCODE><ADAPTER><SEQUENCE><RC_ADAPTER>, where ADAPTER is defined in the adapters file,
        and RC_ADAPTER is defined in the rcadapters file.

    :param input_f: Filepath to input file or folder.
    :param adapters: Filepath to adapters file.
    :param adaptersrc: Filepath to reverse complemented adapters file.
    :param outdir: Filepath to output directory.
    :param threads: Number of processes to use to trim the input fileset.
    :param aux_params: A dictionary of program-specific named-parameters.
    """

    makeDirOrdie(outdir)
    inputs = getInputFiles(input_f)
    pool = init_pool(min(len(inputs), threads))
    debugPrintInputInfo(inputs, "trim adapters from")

    if program == "flexbar":
        clean_trim_adapters_flexbar(inputs, adapters, adaptersrc, outdir, aux_params["allowedns"], pool)
    cleanup_pool(pool)
