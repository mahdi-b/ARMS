import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq


def removeMacseRefs(file_to_clean, reference_fie, output_file_name):
    """Remove the reference sequences from the MACSE files and remove the non nucleotide characters from the sequences.
       we need the database seq. names to remove them from the results files

    :param file_to_clean: Filepath to the list of sequences to clean.
    :param reference_fie: The file containing all reference sequences used to align file_to_clean.
    :param output_file_name: Filepath to where the cleaned sequences should be written.

    :return Filepath to the output file.
    """

    good_seqs = []
    if os.path.isfile(output_file_name):
        os.remove(output_file_name)

    output = open(output_file_name, 'w')
    print "Parsing the DB"
    db_seq_names = SeqIO.to_dict(SeqIO.parse(reference_fie, "fasta")).keys()

    print "Filtering reference sequences"
    i = 0
    for mySeq in SeqIO.parse(file_to_clean, 'fasta'):
        if mySeq.id not in db_seq_names:
            mySeq.seq = Seq(str(mySeq.seq[2:]).replace("-", ""))  # remove the !! from the beginning
            good_seqs.append(mySeq)
            i +=1
            if i % 5000 ==0:
                SeqIO.write(good_seqs, output, 'fasta')
                good_seqs = []
    SeqIO.write(good_seqs, output, 'fasta')
    return output_file_name


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: MACSE_output_file reference_file output_file "
    else:
        removeMacseRefs(*sys.argv[1:4])
