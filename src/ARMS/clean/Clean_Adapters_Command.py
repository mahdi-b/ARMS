from classes.ChewbaccaCommand import *
from Clean_Adapters_Program_Flexbar import Clean_Adapters_Program_Flexbar


class Clean_Adapters_Command(ChewbaccaCommand):
    """Removes sequencing adapters (and preceeding barcodes) from sequences in input file(s).  Sequences should be in
        the following format: <BARCODE><ADAPTER><SEQUENCE><RC_ADAPTER>.  Valid ADAPTER sequences, and their
        reverse-complements (RC_ADAPTER) should be defined separately in a pair of fasta-formatted files.  Sequences
        passed to this command should have already been demultiplexed, as this process will remove the identifying
        barcode sequences.
    """
    supported_programs = [Clean_Adapters_Program_Flexbar]
    default_program = Clean_Adapters_Program_Flexbar
    command_name = "Clean Adapters"
