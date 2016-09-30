from classes.ChewbaccaCommand import *

from Clean_Adapters_Program_Flexbar import Clean_Adapters_Program_Flexbar


class Clean_Adapters_Command(ChewbaccaCommand):
    """Removes adapters (and preceeding barcodes) from sequences in input file(s).
    """
    supported_programs = [Clean_Adapters_Program_Flexbar]
    default_program = Clean_Adapters_Program_Flexbar
    command_name = "Clean Adapters"
