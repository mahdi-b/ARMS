from classes.ChewbaccaCommand import *

from Demux_Program_Fastx import Demux_Program_Fastx


class Demux_Command(ChewbaccaCommand):
    """Given a set of files, each file is assigned a file_ID#.  Each file is then split into separate child files where
        each file holds only sequences belonging to a single sample.  These child files are named using the sample name
        of the sequences it lists, and the file_ID# of the file it came from.  Demuxing is based on the nucleotide
        barcode prefixing each sequence.
    """
    supported_programs = [Demux_Program_Fastx]
    default_program = Demux_Program_Fastx
    command_name = "Demux"
