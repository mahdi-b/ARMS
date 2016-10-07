from classes.ChewbaccaCommand import *

from Dereplicate_Program_Vsearch import Dereplicate_Program_Vsearch


class Dereplicate_Command(ChewbaccaCommand):
    """Dereplicates a fasta file by grouping identical reads together under one representative sequence.  The
    number of duplicate/replicant sequences each representative sequence represents is given by a 'replication count' at
    the end of
    the sequence name in output fasta file.  If a .groups file is provided, then previous replication counts will
    be take in into account (e.g. Imagine a representative sequence X that represents 3 sequences.  If X is found to be
    a replicant of another sequence Y, X will add 3 to replication count of Y).  Replication counts are denoted with a
    suffix of '_K' on the sequence name, where K is the replication count for the group that sequence represents.


    **Inputs**:
        * One or more fasta files to dereplicate.
        * Optional: :ref:`.groups file` - A list of representative names and the names of their replicant \
                                            sequences.  You likely have one of these files if you've previously run a \
                                            clustering or dereplication command.

    **Outputs**:
        * _counts.fasta file - A fasta file with unique sequences and their replication counts.
        * _derep:ref:`.groups file` - A list of representative names and the names of their replicant \
                                            sequences.

    **Notes**:
        * This command only dereplicates within each fasta file (not across all files). \
            This means a sequence in one file will be unique within that file, but might exist in another file. \
            To ensure sequences are uniqe across an entire dataset, merge all fasta files into one file, then \
            dereplicate that fasta file.


        * Each input file will generate a corresponding _count file.
        * If an input .groups file is not provided, then each input fasta file will generate a new groups file named \
            <file_name>_derep.groups.  If an input .groups file IS provided, then a single groups file named \
            'dereplicated_updated.groups' will be generated.
        * The output .groups file is needed by downstream Chewbacca processes (Dereplication, Clustering,\
                Building the OTU Table).
        * The order of sequence names in the \*_counts.fasta and .groups file is arbitrary.

    **Example**:

    ::

        ./
            Data.fasta
                >seq1
                AAA
                >seq2
                AAA
                >seq3
                AAAG
                >seq4
                AAAGT
                >seq7
                AAAGT

            test.groups
                seq4	seq4 seq5 seq6

    In the above example, test.groups indicates that seq4 is a sequence that has previously been identified as a
    representative (in some earlier round of clustering or dereplication).  Note that seq4 is a representative
    for a 'group' of identical sequences and therefore listed within that group.

    ``$ python chewbacca.py dereplicate_fasta -i Data.fasta -o rslt -g test.groups``

    ::

        rslt/Data_counts.fasta:
            >seq4_4
            AAAGT
            >seq1_2
            AAA
            >seq3_1
            AAAG

        rslt_groups_files/*.groups:
            seq3	seq3
            seq1	seq2 seq1
            seq4	seq7 seq6 seq5 seq4

    Notice that Data_counts.fasta lists the unique sequences from Data.fasta, and their replication counts.  Also
    notice that seq4 had previous replication data (stored in the Data.groups file).
    """
    supported_programs = [Dereplicate_Program_Vsearch]
    default_program = Dereplicate_Program_Vsearch
    command_name = "Dereplicate"
