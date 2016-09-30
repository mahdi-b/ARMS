from classes.ChewbaccaCommand import *

from Clean_Quality_Program_Trimmomatic import Clean_Quality_Program_Trimmomatic


class Clean_Quality_Command(ChewbaccaCommand):
    """Removes areas of low quality from reads, keeping the longest contiguous segment."""
    supported_programs = [Clean_Quality_Program_Trimmomatic]
    default_program = Clean_Quality_Program_Trimmomatic
    command_name = "Clean_Quality"
