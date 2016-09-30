from classes.ChewbaccaCommand import *

from Cluster_Program_Crop import Cluster_Program_Crop
from Cluster_Program_Swarm import Cluster_Program_Swarm
from Cluster_Program_Vsearch import Cluster_Program_Vsearch


class Cluster_Command(ChewbaccaCommand):
    """Clusters a set of fasta files.  Identical (or identical spanning) sequences are considered \
        replicants.  (100% match).  NOTE: only Clusters within each fasta file (not across all files).  Merge \
        files before hand if you want to dereplciate across multiple files.
    """
    supported_programs = [Cluster_Program_Crop,
                          Cluster_Program_Swarm,
                          Cluster_Program_Vsearch
                          ]
    default_program = Cluster_Program_Swarm
    command_name = "Cluster"
