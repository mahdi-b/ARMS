import sys
from Bio import SeqIO
from countToDict import parseCountFileToCountDict

#TODO DELETE THIS
# Renames sequences by looking up their count in a global count file.
def renameSequencesWithCount(input_fasta, count_file, outfile):
    seeds = parseCountFileToCountDict(count_file)
    print "Parsed countfile to dict"
    i=0
    mySeqs=[]
    # pull the sequence data for each seed from the input fasta
    for mySeq in SeqIO.parse(input_fasta, 'fasta'):
        i+=1

        if seeds.has_key(mySeq.id):
            mySeq.id = "%s_%i" % (mySeq.id, seeds[mySeq.id])
        else:
            mySeq.id = "%s_%i" % (mySeq.id, 1)
        mySeq.name = ""
        mySeq.description =""
        mySeqs.append(mySeq)
    # write a new fasta with only the seed sequences
    SeqIO.write(mySeqs, open(outfile, 'w'),'fasta')

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fastaA fastaB outfile"
        exit()
        renameSequencesWithCount(sys.argv[1:4])
