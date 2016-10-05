from classes.ChewbaccaCommand import *

from Merge_Program_Chewbacca import Merge_Program_Chewbacca


class Merge_Command(ChewbaccaCommand):
    """Concatenates multiple files into a single file.  Useful for combining the results of a parallel operation, or \
        when preparing for cross-sample derepication."""
    supported_programs = [Merge_Program_Chewbacca]
    default_program = Merge_Program_Chewbacca
    command_name = "Merge"
