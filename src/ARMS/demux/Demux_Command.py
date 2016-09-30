from classes.ChewbaccaCommand import *

from Demux_Program_Fastx import Demux_Program_Fastx


class Demux_Command(ChewbaccaCommand):
    """Splits a set of sequences into separate files by sample based on the barcode at the begining of each sequence."""
    supported_programs = [Demux_Program_Fastx]
    default_program = Demux_Program_Fastx
    command_name = "Demux"
