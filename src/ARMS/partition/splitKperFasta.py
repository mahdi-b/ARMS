import sys
from Bio import SeqIO


def splitK(inputFasta, prefix, nbSeqsPerFile, filetype):
    mySeqs = SeqIO.parse(inputFasta, filetype)
    chunk = 0
    sequences = []

    for mySeq in mySeqs:
        mySeq.seq  = mySeq.seq.ungap(".")
        if len(mySeq.seq) < 200:
            continue
        sequences.append(mySeq)
        if len(sequences) % nbSeqsPerFile == 0:
            SeqIO.write(sequences, open("%s_part_%d.%s" % (str(prefix), chunk, filetype), 'w'), filetype)
            sequences=[]
            chunk+=1
    if sequences:
        SeqIO.write(sequences, open("%s_part_%d.%s" % (str(prefix), chunk, filetype), 'w'), filetype)
    print("Split %s into %d parts." % (inputFasta, (chunk + 1)))

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Usage: input_fasta  output_file_prefix  #seqs_per_file  input_filetype"
    else:
        splitK(*sys.argv[1:5])