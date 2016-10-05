from classes.ChewbaccaCommand import *
from util.Convert_Fastq_Fasta_Program_Chewbacca import Convert_Fastq_Fasta_Program_Chewbacca


class Convert_Fastq_Fasta_Command(ChewbaccaCommand):
    """Converts a Fasta-formatted file to a FastQ-formatted file.  Useful for reducing data size and preparing for
        fasta-only operations."""
    supported_programs = [Convert_Fastq_Fasta_Program_Chewbacca]
    default_program = Convert_Fastq_Fasta_Program_Chewbacca
    command_name = "Convert Fastq to Fasta"
