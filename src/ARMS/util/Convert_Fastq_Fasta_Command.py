from classes.ChewbaccaCommand import ChewbaccaCommand
from util.Convert_Fastq_Fasta_Program_Chewbacca import Convert_Fastq_Fasta_Program_Chewbacca


class Convert_Fastq_Fasta_Command(ChewbaccaCommand):
    """Converts a Fasta-formatted file to a FastQ-formatted file.  Useful for reducing data size and preparing for
        fasta-only operations.

        **Inputs**:
            * One or more fastq files to convert to fasta.

        **Outputs**:
            * <filename>.fasta file(s) - Converted fasta files.

        **Example**:

        ::

            ./
                Data.fasta:
                    @Data_ID#1
                    AGACGCGGWACWGGWTGAACWGTWTAYCCYCCATCGATCGATCGTGRTTYTTYGGNCAYCCNGARGTNTA


        ``$ python chewbacca.py trim_adapters  -i Data.fasta -o rslt -a Data.adapters -arc Data.adapters_RC``


        ::

            rslt/
                Data_debarcoded.fastq:
                    @Data_ID#1
                    ATCGATCGATCG
    """
    supported_programs = [Convert_Fastq_Fasta_Program_Chewbacca]
    default_program = Convert_Fastq_Fasta_Program_Chewbacca
    command_name = "Convert Fastq to Fasta"
