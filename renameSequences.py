from Bio import SeqIO
import sys

# Renames sequences for swarm


seeds ={}
for line in open("ALL_PARSED.out", 'r'):
    line = line.rstrip()
    seeds[line.split("\t")[0]] = line.split("\t")[1]

i=0
mySeqs=[]
for mySeq in SeqIO.parse("MACSEOUT_MERGED.fasta", 'fasta'):
    i+=1

    if seeds.has_key(mySeq.id):
        mySeq.id = "%i_%i" % (i, len(seeds[mySeq.id].split(" ")))
    else:
        mySeq.id = "%i_%i" % (i, 1)
    mySeqs.append(mySeq)

SeqIO.write(mySeqs, open(sys.argv[1], 'w'),'fasta')
    
