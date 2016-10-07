from classes.ChewbaccaCommand import *

from Annotate_OTU_Table_Program_Chewbacca import Annotate_OTU_Table_Program_Chewbacca


class Annotate_OTU_Table_Command(ChewbaccaCommand):
    """Annotates an OTU table with Taxonomic names by replacing sequence names in the OTU table with their identified
    taxonomies.

   **Inputs**:
        * An :ref:`OTU table` to annotate.
        * One or more :ref:`.tax` to read annotations from.

    **Outputs**:
        * An OTU table with sequence names replaced by taxonomic names in the input .tax file.

    **Notes**:
        * The input annotation file(s) should list only one identification per sequence name.  If you find more than \
            one taxonomic identity for a sequence, choose only one to include in the input .tax file(s).

    **Example**:

    ::

        ./
            matrix.txt
                OTU	Sample1	Sample2	Sample3	Sample4
                seq3	1	2	0	0
                seq1	2	1	1	0

            data.tax:
                seq1	94483305	99.4	173	55.4    Chordata:Mammalia:Primates:Hominidae:Homo:Homo sapiens

    ``$ python chewbacca.py annotate_matrix -i matrix.txt -a data.tax -o rslt``

    ::

        rslt/
            matrix.txt
                OTU	Sample1	Sample2	Sample3	Sample4
                seq3	1	2	0	0
                Chordata:Mammalia:Primates:Hominidae:Homo:Homo sapiens	2	1	1	0

    """
    supported_programs = [Annotate_OTU_Table_Program_Chewbacca]
    default_program = Annotate_OTU_Table_Program_Chewbacca
    command_name = "Annotate OTU Table"
