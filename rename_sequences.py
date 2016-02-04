from Bio import SeqIO
import sys
import os

prefix = sys.argv[1].split(".")[0]
fileType = sys.argv[3]
if os.path.exists(sys.argv[2]):
    sys.exit("\n---------------------------------\n\nSequence %s exists!\n---------------------------------\n\n")
outFile = open(sys.argv[2], 'a')
i =0
for s in SeqIO.parse(sys.argv[1], fileType):
    s.id ="%sID%s" % (prefix, i)
    s.description =""
    SeqIO.write(s, outFile, fileType)
    i+=1
outFile.close()

