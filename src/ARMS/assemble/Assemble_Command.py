from classes.ChewbaccaCommand import *

from Assemble_Program_Pear import Assemble_Program_Pear


class Assemble_Command(ChewbaccaCommand):
    """Assembles reads from two (left and right) fastq files/directories.  For a set of k forward read files, and k
        reverse read files, return k assembled files.  Matching forward and reverse files should be identically named,
        except for a <forward>/<reverse> suffix that indicates the read orientation.  Two suffix conventions are
        supported:
        '_forwards'/'_reverse',
        and
        '_R1'/'R2'
        Choose ONE suffix style and stick to it.
        e.g. Sample_100_forwards.fq and Sample_100_reverse.fq will be assembled into Sample_100_assembled.fq.
        Alternatively, Sample_100_R1.fq and Sample_100_R2.fq will be assembled into Sample_100_assembled.fq.
        You can provide as many pairs of files as you wish as long as they follow exactly one of the above naming
        conventions.  If a 'name' parameter is provided, it will be used as a filename (not path) prefix for all
        assembled sequence files.
    """
    supported_programs = [Assemble_Program_Pear]
    default_program = Assemble_Program_Pear
    command_name = "Assemble"
