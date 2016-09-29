import os

from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *

from classes.Helpers import *


class Clean_Adapters_Program_Flexbar(ChewbaccaProgram):
    name = "flexbar"

    def execute_program(self):
        args = self.args
        self.clean_trim_adapters_flexbar(args.input_f, args.adapters, args.adaptersrc, args.outdir, args.allowedns,
                                         args.threads)


    def clean_trim_adapters_flexbar(self, input_f, adapters, adaptersrc, outdir, allowedns, threads):
        """Use flexbar to trim adapters and barcodes from sequences.  By default, Flexbar does not allow any 'N' \
            characters in SEQUENCE, and will toss any sequences that do contain 'N'.  To avoid this, use the -u or \
            --allowedns flags to specify the maximum number of 'N's to allow

        :param input_f: Filepath to input file or folder.
        :param adapters: Filepath to a list of adapters.
        :param adaptersrc: Filepath to a list of reverse-complemented adapters.
        :param outdir: Filepath to the output directory.
        :param allowedns: Non-negative integer value indicating the maximum number of 'N's to tolerate in a sequence.
        """
        makeDirOrdie(outdir)
        inputs = getInputFiles(input_f)
        pool = init_pool(min(len(inputs), threads))
        debugPrintInputInfo(inputs, "trim adapters from")
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