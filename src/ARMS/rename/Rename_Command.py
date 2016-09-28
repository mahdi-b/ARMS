from classes.ChewbaccaCommand import *
from Rename_Program_Chewbacca import Rename_Program_Chewbacca


class Rename_Command(ChewbaccaCommand):

    supported_programs = [Rename_Program_Chewbacca]
    default_program = Rename_Program_Chewbacca
    command_name = "Rename"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
