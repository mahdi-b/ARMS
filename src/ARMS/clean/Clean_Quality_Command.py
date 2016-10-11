from classes.ChewbaccaCommand import ChewbaccaCommand
from Clean_Quality_Program_Trimmomatic import Clean_Quality_Program_Trimmomatic


class Clean_Quality_Command(ChewbaccaCommand):
    """Removes regions of low quality from fastq-formatted reads.  These regions are likely sources of error, and would
        be detrimental to other analytical process.  Input sequences to this command should have already been
        demultiplexed, and had their barcodes/adapters removed.  Otherwise, the partial removal of these markers would
        leave behind invalid fragments that would be difficult to detect and likely cause errors down-stream.


        **Inputs**:
            * One or more fastq files to clean.

        **Outputs**:
            * <filename>_cleaned.fastq file(s) - Fastq files, containing sequences with areas of low quality removed.

        **Notes**:
            * Be aware of the program-specific details around 'N' nucleotide characters.
            * Be aware of the program-specific defaults for minimum surviving sequence lengths.

        **Example**:

        ::

            ./
                Data.fasta:
                    @Data_ID#1
                    GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCTTTACAG
                    +
                    !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz%%%zzzz

        The command below asks Chewbacca to trim away any section of length 3 NT in Data_ID#1 that has quality worse \
            than 20, keeping the \
            longer of the remaining ends.  If the remaining sequence at the end of this process is shorter than 15 NT, \
            discard the whole sequence (these values are chosen for illustrative purposes).

        ``$ python chewbacca.py clean_seqs -i Data.fasta -o rslt -m 15 -w 3 -q 20``

        Note that the 'TTT' subsequence has been cut, because its average quality (5) is less than the threshold (20). \
        After this cut, the longest surviving subsequence (the subsequence to the left of the cut) was kept, and the \
        shorter subsequence (to the right of the cut) was discarded.  Because the final sequence is longer than 15NT, \
        it is kept and written to the output file.

        ::

            rslt/
                Data_cleaned.fastq:
                    @Data_ID#1
                    GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTC
                    +
                    !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
    """
    supported_programs = [Clean_Quality_Program_Trimmomatic]
    default_program = Clean_Quality_Program_Trimmomatic
    command_name = "Clean_Quality"
