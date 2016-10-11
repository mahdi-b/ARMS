from classes.ChewbaccaCommand import ChewbaccaCommand
from Assemble_Program_Pear import Assemble_Program_Pear


class Assemble_Command(ChewbaccaCommand):
    """Assembles reads from two (forward and reverse) fastq files/directories.  For a set of k forward read files, and k
        reverse read files, return k assembled files.  Matching forward and reverse files should be identically named,
        except for a <forward>/<reverse> suffix that indicates the read orientation.  The two suffix conventions below \
        are supported.  Choose ONE suffix style and stick to it!  Mixed suffixes are not supported.

        ::

            _forwards/_reverse
            and
            _R1/_R2


        **Inputs**:
            * fastq file(s) with left reads
            * fastq file(s) with right reads

        **Outputs**:
            * fastq File(s) with assembled reads

        **Notes**:
            * Choose ONE suffix style and stick to it!  Mixed suffixes are not supported. \
                e.g.\
                Sample_100_forwards.fq and Sample_100_reverse.fq will be assembled into Sample_100_assembled.fq. \
                Simmilarly, Sample_100_R1.fq and Sample_100_R2.fq will be assembled into Sample_100_assembled.fq. \
                However, Sample_100_forwards.fq and Sample_100_R2.fq are not guaranteed to be matched.

            * You can provide as many pairs of files as you wish as long as they follow exactly one of the above  \
                naming \
                conventions.  If a 'name' parameter is provided, it will be used as a filename (not path) prefix for \
                 all \
                assembled sequence files.

        **Example**

        Assuming a forwards read file 'Data_R1.fq' and a reverse reads file 'Data_R1.fq',

        ::

            ./
                Data_R1.fq
                Data_R2.fq

        ``$ python chewbacca.py assemble -n BALI -f Data_R1.fq  -r Data_R2.fq  -o rslt``

        ::

            rslt/
                BALI_DATA.assembled.fq
    """
    supported_programs = [Assemble_Program_Pear]
    default_program = Assemble_Program_Pear
    command_name = "Assemble"
