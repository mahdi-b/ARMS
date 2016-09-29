from classes.ChewbaccaCommand import *
from otu.Build_OTU_Table_Program_Chewbacca import Build_OTU_Table_Program_Chewbacca


class Build_OTU_Table_Command(ChewbaccaCommand):

    supported_programs = [Build_OTU_Table_Program_Chewbacca]
    default_program = Build_OTU_Table_Program_Chewbacca
    command_name = "Build OTU Table"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
