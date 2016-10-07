from classes.ChewbaccaCommand import *
from otu.Build_OTU_Table_Program_Chewbacca import Build_OTU_Table_Program_Chewbacca


class Build_OTU_Table_Command(ChewbaccaCommand):
    """Builds an OTU table using a .groups, .samples, and .barcodes file.  The OTU table shows OTU (group) abundance by
    sample.

    **Inputs**:
        * One or more :ref:`.samples`.
        * One or more :ref:`.barcodes`.
        * one or more :ref:`.groups`.

    **Outputs**:
        * matrix.txt - A tab-delimited table mapping OTUs (groups) to their abundance in each sample.

    **Notes**:
        * A sequence name may not appear in more than one group file (or more than one line in a gropus file for \
            that matter!).

    **Example**:

    ::

        ./
            test.barcodes
                Sample1	aaaaaa
                Sample2	aaaaat
                Sample3	aaaaac
                Sample4	aaaaag

            test.groups
                seq3	seq3 seq5 seq4
                seq1	seq2 seq1 seq7 seq6

            test.samples
                seq1	Sample1
                seq2	Sample1
                seq3	Sample1
                seq4	Sample2
                seq5	Sample2
                seq6	Sample2
                seq7	Sample3

    ``$ python chewbacca.py build_matrix -b test.barcodes -g  test.groups -s test.samples -o rslt/``

    ::

        rslt/
            matrix.txt
                OTU	Sample1	Sample2	Sample3	Sample4
                seq3	1	2	0	0
                seq1	2	1	1	0
    """
    supported_programs = [Build_OTU_Table_Program_Chewbacca]
    default_program = Build_OTU_Table_Program_Chewbacca
    command_name = "Build OTU Table"
