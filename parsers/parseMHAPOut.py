import sys
from Bio import SeqIO

def findUnmatchedSeqs(fasta_to_clean, mhap_outfile, unmatched_fasta):
    """Takes in a Mhap file,

    :param fastaA:  The first fasta file.  Preferably the larger fasta file.
    :param fastaB:  The second fasta file. Preferably the smaller fasta file.
    :param outFile: Where to write the resuling fasta to.
    :return:
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
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fasta_to_clean  mhap_outfile  unmatched_fasta"
    else:
        findUnmatchedSeqs(*sys.argv[1:4])
