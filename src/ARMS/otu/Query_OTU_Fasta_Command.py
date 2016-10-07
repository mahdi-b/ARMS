from classes.ChewbaccaCommand import *
from Query_OTU_Fasta_Program_Vsearch import Query_OTU_Fasta_Program_Vsearch


class Query_OTU_Fasta_Command(ChewbaccaCommand):
    """Aligns sequences in a fasta file against those in a reference fasta in order to determine OTU identity.

    **Inputs**:
        * One or more fasta files containing sequences to identify.
        * A curated fasta file of high quality sequences and known species.
        * A two-column, tab-delimited text file mapping sequence names in the curated fasta file to taxonomic \
            identifiers.

    **Outputs**:
        * A `.tax file`- A five-column, tab-delimited output file listing the *'query', 'target', 'id', 'alnlen'*, and \
            *'qcov'* fields from a blast6 output file. Maps input sequences to their identified Taxonomic names.

    **Notes**:
        * The files 'bold.fna' and 'seq_lin.mapping' are included in the Chewbacca Docker distributions.

    **Example**:

    ::

        ~/ARMS/data/
            bold.fna # A precompiled fasta file of data from BOLD.
                >GBMAA1117-14
                GGGCTTTTGCGGGTATGATAGGAACAGCATTTAGTATGCTTATTAGGTTAGAACTATCTTCCCCAGGGTCTATGTTAGGAGATGATCATTTATATAAT
                GTTATAGTAACAGCTCATGCATTTGTAATGATATTTTTTTTAGTTATGCCAGTAATGATTGGGGGTTTTGGTAATTGGTTAGTACCTTTATATATTGG
                TGCCCCGGATATGGCTTTTCCTAGATTAAATAATATTAGTTTTTGGTTATTACCTCCGGCGCTTACTTTATTATTAGGTTCGGCTTTTGTAGAACAAG
                GGGCTGGGACAGGTTGGACAGTTTATCCGCCTTTATTTAGTATTCAAACTCATTCTGGGGGGTCTGTGGATATGGTAATATTTAGTTTACATTTAGCT
                GGAATATCTTCTATATTAGGGGCTATGAATTTTATAACAACAATCTTTAATATGAGGTCTCCGGGAGTAACTATGGATAGAATGCCTTTATTTGTTTG
                ATCTGTTTTAGTAACTGCTTTTTTATTATTATTATCATTGCCAGTATTAGCTGGTGCCATAACAAGTCTTTTAACCGATCGAGATTTTAATACTACAT
                TT

            seq_lin.mapping # A precompiled two-column tab file of (Taxa) for the entries in 'bold.fna'.
                GBMAA1117-14	Animalia;Porifera;Demospongiae;Haplosclerida;Phloeodictyidae;;Calyx;Calyx podatypa

        ./
            Data.fasta:
                >seq1
                ACTATCAGGCATTCAAGCCCATTCAGGGGGAGCAGTAGATATGGCTATATTTAGTCTACATCTAGCTGGTGTATCCTCTATTTTAAGTTCTATAAACT
                TTATAACTACTATAATTAATATGAGGGTTCCTGGGATGAGTATGCATAGATTACCTCTATTCGTATGGTCTGTATTAGTTACTACAATATTATTGTTG
                TTATCTTTACCAGTATTAGCTGGTGGAATTACAATGTTATTGACAGATAGAAATTTTAATACAACATTCTTTGACCCTGCGGGAGGAGGAGATCCTAT
                TTTATTCCAGCACTTATTT

    ``$ python chewbacca.py query_fasta -i Data.fasta -o rslt -r ~/ARMS/data/bold.fna -x ~/ARMS/data/seq_lin.mapping``

    ::

        rslt/
                Data_result.out
                    seq1	GBMAA1117-14	90.6	265	84.7 Animalia;Porifera;Demospongiae;Haplosclerida;Phloeodictyidae;;Calyx;Calyx podatypa
    """
    supported_programs = [Query_OTU_Fasta_Program_Vsearch]
    default_program = Query_OTU_Fasta_Program_Vsearch
    command_name = "Query OTU Identity"


    def execute_command(self):
        self.get_program(self.args.program).execute_program()
