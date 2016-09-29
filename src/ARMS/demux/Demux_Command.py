from classes.ChewbaccaCommand import *

from Demux_Program_Fastx import Demux_Program_Fastx


class Demux_Command(ChewbaccaCommand):

    supported_programs = [Demux_Program_Fastx]
    default_program = Demux_Program_Fastx
    command_name = "Demux"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
