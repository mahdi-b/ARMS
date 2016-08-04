import sys
from Bio import SeqIO

def intersection(fastaA, fastaB, outFile):
    """Writes the intersection of fastaA and fastaB to outFile.
        i.e. the records for sequences in both fastaA and fasta B.

    :param fastaA:  The first fasta file.  Preferably the smaller fasta file.
    :param fastaB:  The second fasta file. Preferably the larger fasta file.
    :param outFile: Where to write the resuling fasta to.
    :return:
    """
    commonRecords = {}
    with open(outFile,'w') as out:
        for record in SeqIO.parse(fastaA, "fasta"):
            commonRecords[record.id]=""

        foundIDs = commonRecords.keys()
        for record in SeqIO.parse(fastaB, "fasta"):
            if record.id in foundIDs:
                out.write(record.format("fasta"))
    out.close()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fastaA fastaB outfile"
        exit()
    intersection(sys.argv[1], sys.argv[2], sys.argv[3])
