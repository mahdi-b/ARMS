from classes.ChewbaccaCommand import ChewbaccaCommand
from otu.Query_OTU_DB_Program_Vsearch import Query_OTU_DB_Program_Vsearch


class Query_OTU_DB_Command(ChewbaccaCommand):
    """Aligns sequences in a fasta file against those in a reference database in order to determine OTU identity.


    **Inputs**:
        * One or more fasta files containing sequences to identify.
        * A curated fasta file of high quality sequences and known species.
        * A database containing taxonomic identifiers for sequences in the curated fasta file.

    **Outputs**:
        * A :ref:`.tax`.

    **Notes**:
        * The files COI.fasta and ncbi.db are included in the Chewbacca Docker distributions.

    **Example**:

    ::

        ~/ARMS/refs/

            COI.fasta # A precompiled fasta file of COI data from NCBI.
                >94483305
                AGGACGGATCAGACGAAGAGGGGCGTTTGGTATTGGGTTATGGCAGGGGGTTTTATATTGATAATTGTTGTGATGAAATT
                GATGGCCCCTAAGATAGAGGAGACACCTGCTAGGTGTAAGGAGAAGATGGTTAGGTCTACGGAGGCTCCAGGGTGGGAGT

            ncbi.db # A precompiled database of (Taxa) for the entries in 'COI.fasta'.


        data/
            Data.fasta:
                >seq1
                GAATAGGTGTTGGTATAGAATGGGGTCTCCTCCTCCGGCGGGGTCGAAGAAGGTGGTGTTGAGGTTGCGGTCTGTTAGTAGTATAGTGATGCCAGCAG
                CTAGGACTGGGAGAGATAGGAGAAGTAGGACTGCTGTGATTAGGACGGATCAGACGAAGAGGGGCGTTTGGTATTGGGTTATGGCAGGGGGTTTTATA
                TTGATAATTGTTGTGAGGAAATTGATGGCCCCTAAGATAGAGGAGACACCTGCTAGGTGTAAGGAGAAGATGGTTAGGTCTACGGAGGCTCCAGGGTG
                GGAGTAGTTCCCTGCTAA

    ``$ python chewbacca.py query_db -i Data.fasta -o out/ -r ~/ARMS/refs/COI.fasta -d ~/ARMS/refs/ncbi.db``

    ::

        rslt/
            Data_result.out
                seq1	94483305	99.4	173	55.4    Chordata:Mammalia:Primates:Hominidae:Homo:Homo sapiens
    """
    supported_programs = [Query_OTU_DB_Program_Vsearch]
    default_program = Query_OTU_DB_Program_Vsearch
    command_name = "Query OTU Identity"
