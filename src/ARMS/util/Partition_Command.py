from classes.ChewbaccaCommand import *
from util.Partition_Program_Chewbacca import Partition_Program_Chewbacca


class Partition_Command(ChewbaccaCommand):
    """A utility command that partitions a fasta/fastq file into  a set of files (of the same file format), with a \
        user-specified (maximum) number of sequences per file.  Allows users to partition a large file into segments, \
        and perform discrete operations in parallel over those segments.
    """
    supported_programs = [Partition_Program_Chewbacca]
    default_program = Partition_Program_Chewbacca
    command_name = "Partition"
