from classes.ChewbaccaCommand import ChewbaccaCommand
from Demux_Name_Program_Chewbacca import Demux_Program_Chewbacca
from Demux_Barcode_Program_Fastx import Demux_Program_Fastx


class Demux_Name_Command(ChewbaccaCommand):
    """Given a set of files, each file is assigned a file_ID#.  Each file is then split into separate child files where
        each file holds only sequences belonging to a single sample.  These child files are named using the sample name
        for the sequences it lists, and the file_ID# of the file it came from.  Demuxing is based on unique
        sample names contained in sequence names.

    **Inputs**:
        * One or more fasta/fastq files to demux.  Sequences in these files should be named with the sample they came
            from.
        * A single .barcodes file: A :ref:`.barcodes`, listing samples as they appear in sequence names, but actual \
            barcode sequences can be made up.  This command will only make use of barcode names.

    **Outputs**:
        * <sample_name>_<file_id#>_ demux.<ext> file(s) - <fasta/fastq> files, containing all the sequences from file \
            <file_id#>, which had a sequence name containing sample <sample_name>.
        * unmatched_<file_id#>_ demux.<ext> file(s) - <fasta/fastq> files, containing sequences from file \
            <file_id#>, whose barcode did not match any of those listed in the .barcodes file.

    **Notes**:
        * The assignment of ID# to file should be treated as an arbitrary process and should not used for record \
            keeping.
        * Each input file will generate its own unmatched_* file (if applicable).

    **Example**:

    ::

        data/
            Data1.fasta:
                @SampleA:001
                AAAAAAAAAAAA
                @SampleAA:002
                AAAAAAAAAAAT
                @SampleA1:003
                AAAAAAAAAAAC
                @Sample_B:001
                AAAAAAAAAAAG

            Data2.fasta:
                @SampleAA:001
                GAAAAAAAAAAA
                @SampleA:002
                TAAAAAAAAAAA
                @Seq8
                CAAAAAAAAAAA
        ./
            Data.barcodes:
                SampleA         AAA
                SampleAA        AAA
                Sample_B        AAA

    ``$ python chewbacca.py demux_names -i data/ -b Data.barcodes -o rslt``

    Here, we see that Data1.fasta was assigned '0' as an ID#, while Data2.fasta was assigned '1' as an ID#.  Because \
    both files had sequences from SampleA, the sequences from Data1.fasta  were written to SampleA_0_demux.fastq, \
    and those sequences from Data2.fasta were written to SampleA_1_demux.fastq.  The same is true for SampleB.

    ::

        rslt/
            SampleA_0_demux.fastq:
                @SampleA:001
                AAAAAAAAAAAA
                @SampleA1:003
                AAAAAAAAAAAC

            SampleAA_0_demux.fastq:
                @SampleAA:002
                AAAAAAAAAAAT

            SampleB_0_demux.fastq:
                @Sample_B:001
                AAAAAAAAAAAG

            SampleA_1_demux.fastq:
                @SampleA:002
                TAAAAAAAAAAA

            SampleAA_1_demux.fastq:
                @SampleAA:001
                GAAAAAAAAAAA

        rslt_aux/
            unmatched_1_demux.fastq:
                @Seq8
                CGTGTAAAAAAG
    """
    supported_programs = [Demux_Program_Chewbacca]
    default_program = Demux_Program_Chewbacca
    command_name = "Demux"
