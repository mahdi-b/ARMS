from Bio import SeqIO
import sys

# Renames sequences for swarm


seeds ={}
for line in open(sys.argv[2], 'r'):
    line = line.rstrip()
    seeds[line.split()[0]] = line.split("\t")[1]

i=0
mySeqs=[]
for mySeq in SeqIO.parse(sys.argv[1], 'fasta'):
    if seeds.has_key(mySeq.id):
        mySeq.id = "%s_%s" % ( mySeq.id, len(seeds[mySeq.id].split(" ")))
        mySeq.description=""
    else:
        mySeq.id = "%s_%s" % (mySeq.id,1)
        mySeq.description=""
    mySeqs.append(mySeq)

SeqIO.write(mySeqs, open(sys.argv[3], 'w'),'fasta')
    
