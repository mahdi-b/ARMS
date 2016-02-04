from Bio import SeqIO
import sys

mySeqs = SeqIO.parse(sys.argv[1], sys.argv[4])
prefix = sys.argv[2]
nbSeqsPerFile = int(sys.argv[3])
chunk = 0
sequences = []

for mySeq in mySeqs:
    mySeq.seq  = mySeq.seq.ungap(".")
    if len(mySeq.seq) < 200:
        continue
    sequences.append(mySeq)
    if len(sequences) % nbSeqsPerFile == 0:
        SeqIO.write(sequences, open(str(prefix)+"_"+str(chunk)+".fa", 'w'), sys.argv[4])
        sequences=[]
        chunk+=1
if sequences:
    SeqIO.write(sequences, open(str(prefix)+str(chunk)+".fa", 'w'), sys.argv[4])
    
        
