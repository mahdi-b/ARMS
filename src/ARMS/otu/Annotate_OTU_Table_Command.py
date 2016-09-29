from Annotate_OTU_Table_Program_Chewbacca import Annotate_OTU_Table_Program_Chewbacca
from classes.ChewbaccaCommand import *


class Annotate_OTU_Table_Command(ChewbaccaCommand):

    supported_programs = [Annotate_OTU_Table_Program_Chewbacca]
    default_program = Annotate_OTU_Table_Program_Chewbacca
    command_name = "Annotate OTU Table"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
