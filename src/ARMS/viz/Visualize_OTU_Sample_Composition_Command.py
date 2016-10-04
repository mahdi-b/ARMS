from classes.ChewbaccaCommand import *
from viz.Visualize_OTU_Sample_Composition_Program_Chewbacca import Visualize_OTU_Sample_Composition_Program_Chewbacca


class Visualize_OTU_Sample_Composition_Command(ChewbaccaCommand):
    supported_programs = [Visualize_OTU_Sample_Composition_Program_Chewbacca]
    default_program = Visualize_OTU_Sample_Composition_Program_Chewbacca
    command_name = "Visualize OTU Sample Composition"
