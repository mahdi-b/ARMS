from align.Align_Clean_Program_Macse import Align_Clean_Program_Macse
from classes.ChewbaccaCommand import *


class Align_Clean_Command(ChewbaccaCommand):
    """Cleans aligned files by removing gap characters and reference sequences from the file.  Sequences passed to this
        command should have previously been aligned.
    """
    supported_programs = [Align_Clean_Program_Macse]
    default_program = Align_Clean_Program_Macse
    command_name = "Align_Clean"
