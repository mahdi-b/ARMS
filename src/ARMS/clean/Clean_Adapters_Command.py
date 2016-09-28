from classes.ChewbaccaCommand import *
from Clean_Adapters_Program_Flexbar import Clean_Adapters_Program_Flexbar


class Clean_Adapters_Command(ChewbaccaCommand):
    """Trim adapters (and preceeding barcodes) from sequences in input file(s).  Sequences should be in
        the following format: <BARCODE><ADAPTER><SEQUENCE><RC_ADAPTER>, where ADAPTER is defined in the adapters file,
        and RC_ADAPTER is defined in the rcadapters file.


    :param adapters: Filepath to adapters file.
    :param adaptersrc: Filepath to reverse complemented adapters file.
    :param outdir: Filepath to output directory.
    :param threads: Number of processes to use to trim the input fileset.
    :param aux_params: A dictionary of program-specific named-parameters.
    """
    supported_programs = [Clean_Adapters_Program_Flexbar]
    default_program = Clean_Adapters_Program_Flexbar
    command_name = "Clean Adapters"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
