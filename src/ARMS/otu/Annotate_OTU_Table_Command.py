from classes.ChewbaccaCommand import *

from Annotate_OTU_Table_Program_Chewbacca import Annotate_OTU_Table_Program_Chewbacca


class Annotate_OTU_Table_Command(ChewbaccaCommand):

    supported_programs = [Annotate_OTU_Table_Program_Chewbacca]
    default_program = Annotate_OTU_Table_Program_Chewbacca
    command_name = "Annotate OTU Table"
