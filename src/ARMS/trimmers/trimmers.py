from src.ARMS.classes.Helpers import *
from src.ARMS.classes.ProgramRunner import *


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
        trim_adapters_flexbar(inputs, adapters, adaptersrc, outdir, aux_params["allowedns"], pool)
    cleanup_pool(pool)


def trim_adapters_flexbar(inputs, adapters, adaptersrc, outdir, allowedns, pool):
    """Use flexbar to trim adapters and barcodes from sequences.  By default, Flexbar does not allow any 'N' \
        characters in SEQUENCE, and will toss any sequences that do contain 'N'.  To avoid this, use the -u or \
        --allowedns flags to specify the maximum number of 'N's to allow

    :param inputs: A list of input files to trim.
    :param adapters: Filepath to a list of adapters.
    :param adaptersrc: Filepath to a list of reverse-complemented adapters.
    :param outdir: Filepath to the output directory.
    :param allowedns: Non-negative integer value indicating the maximum number of 'N's to tolerate in a sequence.
    :param pool: A fully-initalized multiprocessing.Pool object.
    """

    # "flexbar":  "flexbar -r \"%s\" -t \"%s\" -ae \"%s\" -a \"%s\"",
    printVerbose("Trimming barcodes and adapters with flexbar")
    temp_file_name_template = "%s/temp_%s"
    debarcoded_file_name_template = "%s/%s_debarcoded"
    # Trim adapters from the left
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.TRIM_FLEXBAR,
                                                      [input_file,
                                                       temp_file_name_template % (outdir, strip_ixes(input_file)),
                                                       "LEFT", adapters, allowedns],
                                                      {"exists": [input_file, adapters]})
                                        for input_file in inputs], pool)

    temp_files = getInputFiles(outdir, "temp_*")
    debugPrintInputInfo(temp_files, "trim adapters from")

    # Trim the reverse complemented adapters from the right
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.TRIM_FLEXBAR,
                                                      [input_file,
                                                       debarcoded_file_name_template % (outdir,
                                                                                        strip_ixes(input_file)[5:]),
                                                       "RIGHT", adaptersrc,
                                                       allowedns],
                                                      {"exists": [input_file, adaptersrc]})
                                        for input_file in temp_files], pool)
    printVerbose("Done Trimming sequences.")

    # delete temp files
    for file_ in getInputFiles(outdir, "temp_*"):
        os.remove(file_)