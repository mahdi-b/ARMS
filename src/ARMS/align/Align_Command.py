from classes.ChewbaccaCommand import *

from Align_Program_Macse import Align_Program_Macse


class Align_Command(ChewbaccaCommand):
    """Aligns a set of fasta files.  See Align_Program_x documentation for program-specific details.   Input files to
        this command should first be dereplicated.  Doing so will reduce the total number of alignments required, and
        reduce computation time.
    """
    supported_programs = [Align_Program_Macse]
    default_program = Align_Program_Macse
    command_name = "Align"
