from classes.ChewbaccaCommand import *

from Preclean_Program_Bayeshammer import Preclean_Program_Bayeshammer


class Preclean_Command(ChewbaccaCommand):
    """Fixes sequencing errors.
    """
    supported_programs = [Preclean_Program_Bayeshammer]
    default_program = Preclean_Program_Bayeshammer
    command_name = "Preclean"
