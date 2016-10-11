from classes.ChewbaccaCommand import ChewbaccaCommand
from Rename_Program_Chewbacca import Rename_Program_Chewbacca


class Rename_Command(ChewbaccaCommand):
    """
    Renames sequences in a file with their filename and a serial ID#.  Useful for simplifying complex naming systems \
    into \
    human-readable sequence names.  In order to ensure the correct sample names are preserved, it is reccomended that \
    this command be run immediately after the :ref:`Demux Command`.

    **Inputs**:
        * One or more fasta/fastq files to rename.

    **Outputs**:
        * _renamed.<ext> file - A <fasta/fastq> file with the renamed sequences.
        * .samples file - A :ref:`.samples`.
        * .mapping file - A :ref:`.mapping`.

    **Notes**:
        * In order for the .samples file to correctly list the sample name of the sequences in a file, this command \
            should be run immediately after the Demux Command.
        * The *--clip* parameter tells Chewbacca that trailing _<file_ID#> (from the demultiplexing command) should not
            be considered part of the sample name.  By default this is set to True, and should be fine.  If you notice
            parts of your sample names getting clipped off in your .samples file, you should explicitly set this
            parameter to False.
        * Each input file will have a corresponding .samples, .mapping, and _renamed file.
        * The .samples file is needed by downstream Chewbacca processes (Building the OTU Table).
        * The .mapping file is purely for user convenience and record-keeping.

    **Example**:

    ::

        SampleA_0.fasta:
            @M03292:26:000000000-AH6AG:1:1101:22127:1256
            AAAA
            @M03292:26:000000000-AH6AG:1:1101:22127:1257
            AAAT

    ``$ python chewbacca.py rename -i SampleA_0.fasta -o rslt``

    ::

        rslt/SampleA_0_renamed.fasta:
            @SampleA_ID0
            AAAA
            @SampleA_ID1
            AAAT

        rslt_samples/SampleA_0_renamed.samples:
            SampleA_ID0	SampleA
            SampleA_ID1	SampleA

        rslt_aux/SampleA_0_renamed.mapping:
            M03292:26:000000000-AH6AG:1:1101:22127:1256	SampleA_ID0
            M03292:26:000000000-AH6AG:1:1101:22127:1257	SampleA_ID1

    """
    supported_programs = [Rename_Program_Chewbacca]
    default_program = Rename_Program_Chewbacca
    command_name = "Rename"
