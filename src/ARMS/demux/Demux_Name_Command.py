from classes.ChewbaccaCommand import ChewbaccaCommand

from Demux_Barcode_Program_Fastx import Demux_Program_Fastx


class Demux_Command(ChewbaccaCommand):
    """Given a set of files, each file is assigned a file_ID#.  Each file is then split into separate child files where
        each file holds only sequences belonging to a single sample.  These child files are named using the sample name
        for the sequences it lists, and the file_ID# of the file it came from.  Demuxing is based on the nucleotide
        barcode prefixing each sequence.

    **Inputs**:
        * One or more fasta/fastq files to demux.
        * A single .barcodes file: A :ref:`.barcodes`.

    **Outputs**:
        * <sample_name>_<file_id#>_ demux.<ext> file(s) - <fasta/fastq> files, containing all the sequences from file \
            <file_id#>, which had a barcode corresponding to sample <sample_name>.
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
                @Seq4
                AGACGCAAAAAA
                @Seq5
                AGTGTAAAAAAT


            Data2.fasta:
                @Seq6
                AGACGCAAAAAC
                @Seq7
                AGTGTAAAAAAG
                @Seq8
                CGTGTAAAAAAG
        ./
            Data.barcodes:
                SampleA        AGACGC
                SampleB        AGTGTA

    ``$ python chewbacca.py demux_samples -i data/ -b Data.barcodes -o rslt``

    Here, we see that Data1.fasta was assigned '0' as an ID#, while Data2.fasta was assigned '1' as an ID#.  Because \
    both files had sequences from SampleA, the sequences from Data1.fasta  were written to SampleA_0_demux.fastq, \
    and those sequences from Data2.fasta were written to SampleA_1_demux.fastq.  The same is true for SampleB.

    ::

        rslt/
            SampleA_0_demux.fastq:
                @Seq4
                AGACGCAAAAAA

            SampleB_0_demux.fastq:
                @Seq5
                AGTGTAAAAAAT

            SampleA_1_demux.fastq:
                @Seq6
                AGACGCAAAAAC

            SampleB_1_demux.fastq:
                @Seq7
                AGTGTAAAAAAG

        rslt_aux/
            unmatched_0_demux.fastq:
                @Seq8
                CGTGTAAAAAAG
    """
    supported_programs = [Demux_Program_Fastx]
    default_program = Demux_Program_Fastx
    command_name = "Demux"
