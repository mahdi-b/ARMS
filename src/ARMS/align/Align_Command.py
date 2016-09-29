from classes.ChewbaccaCommand import *
from Align_Program_Macse import Align_Program_Macse


class Align_Command(ChewbaccaCommand):
    """Aligns a set of fasta files.  Identical (or identical spanning) sequences are considered \
        replicants.  (100% match).  NOTE: only Aligns within each fasta file (not across all files).  Merge \
        files before hand if you want to dereplciate across multiple files.
    """
    supported_programs = [Align_Program_Macse]
    default_program = Align_Program_Macse
    command_name = "Align"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()


