from classes.ChewbaccaCommand import *

from Dereplicate_Program_Vsearch import Dereplicate_Program_Vsearch


class Dereplicate_Command(ChewbaccaCommand):
    """Dereplicates a set of fasta files.  Identical (or identical spanning) sequences are considered \
        replicants.  (100% match).  NOTE: only dereplicates within each fasta file (not across all files).  Merge \
        files before hand if you want to dereplciate across multiple files.
    """
    supported_programs = [Dereplicate_Program_Vsearch]
    default_program = Dereplicate_Program_Vsearch
    command_name = "Dereplicate"
