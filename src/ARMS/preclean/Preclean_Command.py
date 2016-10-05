from classes.ChewbaccaCommand import *

from Preclean_Program_Bayeshammer import Preclean_Program_Bayeshammer


class Preclean_Command(ChewbaccaCommand):
    """ Attempts to fix minor sequencing errors caused by pyrosequencing.  By reducing errors prior to sequence assembly,
        a greater number of paired reads can be sucessfully assembled.

        **Inputs**:
            * Left reads file(s).
            * Right reads file(s).

        **Outputs**:
            * <Left reads file>_corrected.fastq file(s).
            * <Right reads file>_corrected.fastq file(s).

        **Example**:

        Assuming a forwards read file 'Data_R1.fq' and a reverse reads file 'Data_R1.fq',

        ``$ python chewbacca.py preclean -f Data_R1.fq -r Data_R2.fq  -o rslt``

        ::

                rslt/
                    Data_R1_corrected.fq
                    Data_R1_corrected.fq

        """
    supported_programs = [Preclean_Program_Bayeshammer]
    default_program = Preclean_Program_Bayeshammer
    command_name = "Preclean"
