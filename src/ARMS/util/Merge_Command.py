from classes.ChewbaccaCommand import *

from Merge_Program_Chewbacca import Merge_Program_Chewbacca


class Merge_Command(ChewbaccaCommand):

    supported_programs = [Merge_Program_Chewbacca]
    default_program = Merge_Program_Chewbacca
    command_name = "Merge"