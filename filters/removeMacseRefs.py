import os
from Bio import Seq, SeqIO

# Remove the reference sequences from the MACSE files and remove the non nucleotide characters from the sequences.
#       we need the datbase seq. names to remove them from the results files

def removeMacseRefs(file_to_clean, reference_fie, output_file_name):
    """

    :param file_to_clean: Filepath to the list of sequences to clean.
    :param reference_fie: The file containing all reference sequences used to align file_to_clean.
    :param output_file_name: Filepath to where the cleaned sequences should be written.
    """

    good_seqs = []
    if os.path.isfile(output_file_name):
        os.remove(output_file_name)

    output = open(output_file_name, 'w')
    print "Parsing the DB"
    db_seq_names = SeqIO.to_dict(SeqIO.parse(reference_fie, "fasta")).keys()

    print "Filtering reference seuences"
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