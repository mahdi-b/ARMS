from classes.ChewbaccaCommand import ChewbaccaCommand
from util.Ungap_Program_Chewbacca import Ungap_Program_Chewbacca


class Ungap_Command(ChewbaccaCommand):
    """Removes target characters from a fasta/fastq file.  Useful for removing gap characters from sequence alignments.

    **Inputs**:
        * One or more fasta/fastq files to clean.
        * A string of gap characters to remove.

    **Outputs**:
        * \*_cleaned.<ext> file - A <fasta/fastq> file with gap characters removed from its sequences.

    **Example**:

    ::

        Data.fasta:
            >seq1
            AAAAA.A*A-A

    ``$ python chewbacca.py ungap_fasta -i Data.fasta -o rslt -f fasta -g ".*-"``

    ::

        rslt/Data.fasta:
            >seq1
            AAAAAAAA

    """
    supported_programs = [Ungap_Program_Chewbacca]
    default_program = Ungap_Program_Chewbacca
    command_name = "Ungap"
