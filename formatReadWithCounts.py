from Bio import SeqIO
import operator
import sys
import os

seedSizes = {}

if os.path.isfile(sys.argv[3]):
    os.remove(sys.argv[3])


print "reading uc_parsed.out file"
nbLines = 0
for line in open(sys.argv[2]):
    #line = line.rstrip()
    data = line.rstrip().split()
    size = 0

    for d in data:
        size += int(d.split("_")[1])
    seedSizes[data[0]] = size
    if nbLines % 1000000 == 0:
        print "%s lines processed" % nbLines 
    nbLines +=1

print "Done reading uc_parsed.out file"


seeds = []

print "\nIndexing reads"
reads = SeqIO.index(sys.argv[1], "fasta")
print "Done indexing reads"

print "\nRenaming sequences"
for item in sorted(seedSizes.items(), key=operator.itemgetter(1), reverse=True):
    s = reads[item[0]]

    s.id = "%s_%s" % (item[0].split("_")[0], seedSizes[item[0]])
    s.description = ""
    seeds.append(s)
    if len(seeds) == 500000:
        SeqIO.write(seeds, open(sys.argv[3], 'a'), 'fasta')
        seeds =[]
print "Done renaming sequences"
print "Completed Susccessfully"
