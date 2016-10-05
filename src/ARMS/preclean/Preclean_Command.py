from classes.ChewbaccaCommand import *

from Preclean_Program_Bayeshammer import Preclean_Program_Bayeshammer


class Preclean_Command(ChewbaccaCommand):
    """Attempts to fix minor sequencing errors caused by pyrosequencing.  By reducing errors prior to sequence assembly,
        a greater number of paired reads can be sucessfully assembled.
    """
    supported_programs = [Preclean_Program_Bayeshammer]
    default_program = Preclean_Program_Bayeshammer
    command_name = "Preclean"
