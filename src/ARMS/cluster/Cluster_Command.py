from classes.ChewbaccaCommand import *

from Cluster_Program_Crop import Cluster_Program_Crop
from Cluster_Program_Swarm import Cluster_Program_Swarm
from Cluster_Program_Vsearch import Cluster_Program_Vsearch


class Cluster_Command(ChewbaccaCommand):
    """Clusters a set of fasta files.  This command generates a fasta file of unique sequences
    (each representing a cluster) and a .groups file.  This command also takes an optional .groups file containing
    replication data from previous commands.  If a .groups file is supplied, only one output .groups file is generated
    (regardless of the number of inputs).

    **Inputs**:
        * One or more fasta files to cluster.
        * Optional: :ref:`.groups file` - A list of representative names and the names of their replicant \
                                            sequences.  You likely have one of these files if you've previously run a \
                                            clustering or dereplication command.

    **Outputs**:
        * \*.fasta file - A fasta file with unique sequences and their replication counts.
        * _derep:ref:`.groups file` - A list of representative names and the names of their replicant \
                                            sequences.

    **Notes**:
        * The input fasta file(s) should have been dereplicated before clustering.
        * For a single experiment with multiple fasta files, it is best to merge all input fasta files, dereplicate
            them, then cluster the single merged and dereplicated fasta file.  This provides the best OTU groupings.

    **Example**:

    ::

        ./
            Data.fasta:
                >seq1_3
                AAAAAAAAAA
                >seq2_1
                ATAAAAAAAA
                >seq3_1
                TTTTTTTTTT
                >seq4_1
                TTTTTTATTT
                >seq5_1
                TTTTTTATCT


            Data.groups:
                seq1	seq6 seq1 seq7

    ``$ python chewbacca.py cluster_seqs -i Data.fasta -o rslt -g Data.groups``

    ::

        rslt/
            Data_clustered_seeds.fasta:
                >seq1_4
                AAAAAAAAAA
                >seq3_3
                TTTTTTTTTT

        rslt_groups_files/
            postcluster_updated.groups:
                seq3	seq3 seq5 seq4
                seq1	seq2 seq1 seq7 seq6
    """
    supported_programs = [Cluster_Program_Crop,
                          Cluster_Program_Swarm,
                          Cluster_Program_Vsearch
                          ]
    default_program = Cluster_Program_Swarm
    command_name = "Cluster"
