from classes.ChewbaccaCommand import *
from util.Partition_Program_Chewbacca import Partition_Program_Chewbacca


class Partition_Command(ChewbaccaCommand):
    """A utility command that partitions a fasta/fastq file into  a set of files (of the same file format), with a \
        user-specified (maximum) number of sequences per file.  Allows users to partition a large file into segments, \
        and perform discrete operations in parallel over those segments.

    **Inputs**:
        * One or more fasta/fastq files to partition.
        * C: An integer defining the maximum number of sequences per file

    **Outputs**:
        * <filename>_part_<part_#>.<ext> file(s) - <fasta/fastq> files, with at most C sequences per file.

    **Example**:

    ::

        ./
            Data.fq:
                @Data_ID1
                GATTTGGGG
                +
                !zzzzzzzzz
                @Data_ID2
                GATTTGGGG
                +
                !zzzzzzzzz
                @Data_ID3
                GATTTGGGG
                +
                !zzzzzzzzz


    ``$ python chewbacca.py convert_fastq_to_fasta -i Data.fq -o rslt/``

    ::

        rslt/
            Data.fasta:
                @Data_ID1
                GATTTGGGG
                @Data_ID2
                GATTTGGGG
                @Data_ID3
                GATTTGGGG
    """
    supported_programs = [Partition_Program_Chewbacca]
    default_program = Partition_Program_Chewbacca
    command_name = "Partition"
