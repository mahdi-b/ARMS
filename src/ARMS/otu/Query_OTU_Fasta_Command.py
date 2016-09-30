from classes.ChewbaccaCommand import *

from Query_OTU_Fasta_Program_Vsearch import Query_OTU_Fasta_Program_Vsearch


class Query_OTU_Fasta_Command(ChewbaccaCommand):

    supported_programs = [Query_OTU_Fasta_Program_Vsearch]
    default_program = Query_OTU_Fasta_Program_Vsearch
    command_name = "Query OTU Identity"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
