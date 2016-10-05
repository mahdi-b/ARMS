from classes.ChewbaccaCommand import *

from Dereplicate_Program_Vsearch import Dereplicate_Program_Vsearch


class Dereplicate_Command(ChewbaccaCommand):
    """Dereplicates a set of fasta files.  Identical or spanning (e.g. 'ATA' spans the sequence 'AT') sequences are
        considered replicants.  NOTE: This command only dereplicates within each fasta file (not across all files).
        This means a sequence in one file will be unique within that file, but might exist in another file.
        To ensure sequences are uniqe across an entire dataset, merge all fasta files into one file, then dereplicate
        that fasta file.
    """
    supported_programs = [Dereplicate_Program_Vsearch]
    default_program = Dereplicate_Program_Vsearch
    command_name = "Dereplicate"
