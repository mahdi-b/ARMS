from classes.ChewbaccaCommand import ChewbaccaCommand

from Merge_Program_Chewbacca import Merge_Program_Chewbacca


class Merge_Command(ChewbaccaCommand):
    """Concatenates multiple files into a single file.  Useful for combining the results of a run_parallel operation, or \
        when preparing for cross-sample derepication.

    **Inputs**:
        * A set of files to merge.
        * An <output_filename>.
        * An <output_prefix>.

    **Outputs**:
        * <output_filename>.<output_prefix> - A file consisting of all the input files concatenated together.

    **Notes**:
        * The order of the content in the concatenated files is not guaranteed.

    **Example**:

    ::

        targets/
            Data.fq:
                @Data_ID1
                GATTTGGGG
                +
                !zzzzzzzzz

            Data2.fa:
                @Data_ID1
                GATTTGGGG

            Blah.txt
                Hello World!


    ``$ python chewbacca.py merge_files -i targets/ -o rslt/ -f txt -n Merged``

    ::

        rslt/
            Merged.txt:
                Hello World!
                @Data_ID1
                GATTTGGGG
                +
                !zzzzzzzzz
                @Data_ID1
                GATTTGGGG
    """
    supported_programs = [Merge_Program_Chewbacca]
    default_program = Merge_Program_Chewbacca
    command_name = "Merge"
