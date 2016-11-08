from Bio import SeqIO


def remove_refs_from_macse_out(input_file, db, outfile):
    """Removes reference sequences and gap characters from a macse alignment.

    :param input_file: Filepath to the file to clean.
    :param db: The reference fasta file used in the macse alignment.
    :param outfile: Filepath to the output file.
    """
    dbSeqNames = SeqIO.to_dict(SeqIO.parse(db, "fasta"))
    good_seqs = []
    for mySeq in SeqIO.parse(input_file, 'fasta'):
        # Keep only the non-reference sequences
        if not dbSeqNames.has_key(mySeq.id):
            mySeq.seq = mySeq.seq.ungap("-")[2:]
            good_seqs.append(mySeq)
    SeqIO.write(good_seqs, open(outfile, 'w'), 'fasta')
