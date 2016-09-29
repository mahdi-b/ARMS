from classes.ChewbaccaCommand import *
from align.Align_Clean_Program_Macse import Align_Clean_Program_Macse


class Align_Clean_Command(ChewbaccaCommand):
    """Align_Cleans a set of fasta files.  Identical (or identical spanning) sequences are considered \
        replicants.  (100% match).  NOTE: only Align_Cleans within each fasta file (not across all files).  Merge \
        files before hand if you want to dereplciate across multiple files.
    """
    supported_programs = [Align_Clean_Program_Macse]
    default_program = Align_Clean_Program_Macse
    command_name = "Align_Clean"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()


