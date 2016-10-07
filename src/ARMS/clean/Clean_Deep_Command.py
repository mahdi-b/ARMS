from classes.ChewbaccaCommand import *
from Clean_Deep_Program_Macse import Clean_Deep_Program_Macse


class Clean_Deep_Command(ChewbaccaCommand):
    """Performs an intensive deep-cleaning of sequences to eliminate frameshifts, detect chimeras,
    and determine sequence orientation.  Input files to this command should first be dereplicated.  Doing so will
    reduce the total number of alignments required, and reduce computation time.

    **Inputs**:
        * One or more fasta/fastq files to deep clean.
        * Nucleotide reference fasta.

    **Outputs**:
        * \*_AA - Amino Acid Alignment file, including reference sequences.
        * \*_log.csv - A log listing each input sequence, and deep cleaning results for each sequence.
        * \*_NT -  Nucleotide Alignment file, including reference sequences.

    **Notes**:
        * Sequences that do not meet quality cleaning standards are dropped.
        * The output files contain reference sequences, and odd alignment characters.  Both of these need to be \
            removed by running the Clean_Deep_Repair Command.

    **Example**:

    ::

        Data.fasta
        BIOCODE.fa

    ``$ python chewbacca.py macseAlign -i Data.fasta -o rslt -d BIOCODE.fa``

    ::

        rslt/Data_AA
        rslt/Data_NT
        rslt/Data_log.csv

    """
    supported_programs = [Clean_Deep_Program_Macse]
    default_program = Clean_Deep_Program_Macse
    command_name = "Clean_Deep"
