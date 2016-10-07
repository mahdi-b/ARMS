from classes.ChewbaccaCommand import *
from Clean_Adapters_Program_Flexbar import Clean_Adapters_Program_Flexbar


class Clean_Adapters_Command(ChewbaccaCommand):
    """Removes sequencing adapters (and preceeding barcodes) from sequences in input file(s).  Sequences should be in
        the following format:

        ::

            <BARCODE><ADAPTER><SEQUENCE><RC_ADAPTER>.

        Valid ADAPTER sequences, and their
        reverse-complements (RC_ADAPTER) should be defined separately in a pair of fasta-formatted files.  Sequences
        passed to this command should have already been demultiplexed, as this process will remove the identifying
        barcode sequences.

        **Inputs**:
            * One or more fasta/fastq files to clean.
            * A single :ref:`.adapters` file
            * A single :ref:`.adaptersRC` file

        **Outputs**:
            * <filename>_debarcoded.<ext> file(s) - <fasta/fastq> files, containing sequences with their leading \
                adapters, trailing adapters, and barcodes removed.

        **Notes**:
            * Be aware of the program-specific details around 'N' nucleotide characters.

        **Example**:

        Given Data_ID#1 with  barcode=AGACGC:

        ::

            ./
                Data.fasta:
                    @Data_ID#1
                    AGACGCGGWACWGGWTGAACWGTWTAYCCYCCATCGATCGATCGTGRTTYTTYGGNCAYCCNGARGTNTA

                Data.adapters:
                    >1
                    GGWACWGGWTGAACWGTWTAYCCYCC

                Data.adaptersRC:
                    >first
                    TGRTTYTTYGGNCAYCCNGARGTNTA

        ``$ python chewbacca.py trim_adapters  -i Data.fasta -o rslt -a Data.adapters -arc Data.adapters_RC``


        ::

            rslt/
                Data_debarcoded.fastq:
                    @Data_ID#1
                    ATCGATCGATCG
        """


    supported_programs = [Clean_Adapters_Program_Flexbar]
    default_program = Clean_Adapters_Program_Flexbar
    command_name = "Clean Adapters"
