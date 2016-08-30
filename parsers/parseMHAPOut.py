import sys
from Bio import SeqIO

def findUnmatchedSeqs(fasta_to_clean, mhap_outfile, unmatched_fasta):
    """Given a fasta an the output of running it through mahp (to find matches), remove the sequences that returned
        hits, and write the unmatched sequences to a separate output fasta file.

    :param fasta_to_clean: Filepath to the input fasta that was fed into mhap to generate the out file.
    :param mhap_outfile: Filepath to the .out mhap file generated from running mhap over the input fasta.
    :param unmatched_fasta: Filepath to where the output fasta should be written.
    :return: Filepath to the output unmatched fasta file.
    """
    mhap_rslts = {}
    # parse mhap to dict
    with open(mhap_outfile, 'r') as mhap_file:
        for line in mhap_file:
            data = line.split()
            id = data[0]
            match_data = "\t".join(data[1:])

            if id not in mhap_rslts:
                mhap_rslts[id] = match_data

    # iterate over fasta, if no match, write to outfile
    with open(unmatched_fasta,'w') as out:
        for record in SeqIO.parse(open(fasta_to_clean,'r'), "fasta"):
            if record.id not in mhap_rslts:
                out.write(record.format("fasta"))
    out.close()
    return unmatched_fasta


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fasta_to_clean  mhap_outfile  unmatched_fasta"
    else:
        findUnmatchedSeqs(*sys.argv[1:4])
