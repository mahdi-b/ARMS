from classes.ChewbaccaCommand import *
from util.Convert_Fastq_Fasta_Program_Chewbacca import Convert_Fastq_Fasta_Program_Chewbacca


class Convert_Fastq_Fasta_Command(ChewbaccaCommand):

    supported_programs = [Convert_Fastq_Fasta_Program_Chewbacca]
    default_program = Convert_Fastq_Fasta_Program_Chewbacca
    command_name = "Convert Fastq to Fasta"
