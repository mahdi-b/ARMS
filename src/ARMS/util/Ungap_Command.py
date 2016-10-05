from classes.ChewbaccaCommand import *
from util.Ungap_Program_Chewbacca import Ungap_Program_Chewbacca


class Ungap_Command(ChewbaccaCommand):
    """Removes target characters from a fasta/fastq file.  Useful for removing gap characters from sequence alignments.
    """
    supported_programs = [Ungap_Program_Chewbacca]
    default_program = Ungap_Program_Chewbacca
    command_name = "Ungap"
