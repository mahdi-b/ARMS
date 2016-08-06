import sys
from Bio import SeqIO

def diff(fastaA, fastaB, outFile):
    """Writes the difference of fastaA and fastaB to outFile.  Order matters.
        i.e. Write to file the records in fastaA but not fastaB.

    :param fastaA:  The first fasta file.  Preferably the larger fasta file.
    :param fastaB:  The second fasta file. Preferably the smaller fasta file.
    :param outFile: Where to write the resuling fasta to.
    :return:
    """
    commonRecords = {}
    with open(outFile,'w') as out:
        for record in SeqIO.parse(open(fastaB,'r'), "fasta"):
            commonRecords[record.id]=""
        foundIDs = commonRecords.keys()
        for record in SeqIO.parse(open(fastaA,'r'), "fasta"):
            if not record.id in foundIDs:
                out.write(record.format("fasta"))

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fastaA fastaB outfile"
        exit()
    diff(sys.argv[1], sys.argv[2], sys.argv[3])
