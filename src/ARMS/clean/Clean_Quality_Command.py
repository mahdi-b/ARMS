from classes.ChewbaccaCommand import *

from Clean_Quality_Program_Trimmomatic import Clean_Quality_Program_Trimmomatic


class Clean_Quality_Command(ChewbaccaCommand):

    supported_programs = [Clean_Quality_Program_Trimmomatic]
    default_program = Clean_Quality_Program_Trimmomatic
    command_name = "Clean_Quality"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
