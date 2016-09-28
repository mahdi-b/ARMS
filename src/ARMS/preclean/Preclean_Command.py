from classes.ChewbaccaCommand import *
from Preclean_Program_Bayeshammer import Preclean_Program_Bayeshammer


class Preclean_Command(ChewbaccaCommand):

    supported_programs = [Preclean_Program_Bayeshammer]
    default_program = Preclean_Program_Bayeshammer
    command_name = "Preclean"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
