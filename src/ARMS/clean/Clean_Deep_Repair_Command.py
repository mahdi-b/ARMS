from classes.ChewbaccaCommand import *
from clean.Clean_Deep_Repair_Program_Macse import Clean_Deep_Repair_Program_Macse


class Clean_Deep_Repair_Command(ChewbaccaCommand):
    """Cleans aligned files by removing gap characters and reference sequences from the file.  Sequences passed to this
        command should have previously been aligned.

    **Inputs**:
        * \*_AA - Amino Acid Alignment file, including reference sequences.
        * \*_log.csv - A log listing each input sequence, and deep cleaning results for each sequence.
        * \*_NT -  Nucleotide Alignment file, including reference sequences.
        * Nucleotide reference fasta.
        * \* The original fasta files that were passed in to the Clean_Deep Command
        * \* The Nucleotide reference fasta that was passed to the Clean_Deep Command

    **Outputs**:
        * \*_MERGED.fasta - A clean fasta file with all the surviving sequences from deep cleaning.

    **Notes**:
        * A single *_MERGED.fasta is generated regardless of the number of input files.

    **Example**:

    ::

        BIOCODE.fa

        originalData/Data.fasta

        input/
            Data_AA
            Data_NT
            Data_log.csv

    ``$ python chewbacca.py -i input/ -o out/ -d  BIOCODE.fa -s originalData/``

    ::

        out/
            MACSE_OUT_MERGED.fasta


    """
    supported_programs = [Clean_Deep_Repair_Program_Macse]
    default_program = Clean_Deep_Repair_Program_Macse
    command_name = "Clean_Deep"
