from classes.ChewbaccaCommand import *
from util.Partition_Program_Chewbacca import Partition_Program_Chewbacca


class Partition_Command(ChewbaccaCommand):

    supported_programs = [Partition_Program_Chewbacca]
    default_program = Partition_Program_Chewbacca
    command_name = "Partition"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
