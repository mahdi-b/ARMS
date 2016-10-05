from classes.ChewbaccaCommand import *

from Clean_Quality_Program_Trimmomatic import Clean_Quality_Program_Trimmomatic


class Clean_Quality_Command(ChewbaccaCommand):
    """Removes regions of low quality from fastq-formatted reads.  These regions are likely sources of error, and would
        be detrimental to other analytical process.  Input sequences to this command should have already been
        demultiplexed, and had their barcodes/adapters removed.  Otherwise, partial removal of these markers would be
        difficult to detect and would likely cause errors down-stream.
    """
    supported_programs = [Clean_Quality_Program_Trimmomatic]
    default_program = Clean_Quality_Program_Trimmomatic
    command_name = "Clean_Quality"
