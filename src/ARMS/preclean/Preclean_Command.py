from classes.ChewbaccaCommand import *

from Preclean_Program_Bayeshammer import Preclean_Program_Bayeshammer


class Preclean_Command(ChewbaccaCommand):
    """ Attempts to fix minor sequencing errors caused by pyrosequencing.  By reducing errors prior to sequence assembly,
        a greater number of paired reads can be sucessfully assembled.  Matching forward and reverse files should be \
        identically named,
        except for a <forward>/<reverse> suffix that indicates the read orientation.  The two suffix conventions below \
        are supported. Choose ONE suffix style and stick to it!  Mixed suffixes are not supported.


        ::

            _forwards/_reverse
            and
            _R1/_R2



        **Inputs**:
            * fastq file(s) with left reads .
            * fastq file(s) with right reads .

        **Outputs**:
            * <left reads file>_corrected.fastq file(s).
            * <right reads file>_corrected.fastq file(s).

        **Example**:

        Assuming a forwards read file 'Data_R1.fq' and a reverse reads file 'Data_R1.fq',

        ::

            ./
                Data_R1.fq
                Data_R2.fq

        ``$ python chewbacca.py preclean -f Data_R1.fq -r Data_R2.fq  -o rslt``

        ::

                rslt/
                    Data_R1_corrected.fq
                    Data_R2_corrected.fq

        """
    supported_programs = [Preclean_Program_Bayeshammer]
    default_program = Preclean_Program_Bayeshammer
    command_name = "Preclean"
