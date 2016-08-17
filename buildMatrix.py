import sys
import glob
import os
from collections import defaultdict

seqTosample=dict()
samples =set()

path = os.path.join(sys.argv[1],"*uc_parsed.out")


for f in glob.glob(path):
    sampleName = "_".join(os.path.basename(f).split("_")[1:-4])
    samples.add(sampleName)
    for line in open(f).readlines():
        line = line.rstrip()
        data = line.split("\t")
        seqTosample[data[0]] = sampleName
        if len(data) == 1:
            continue
        for s in data[1].split():
            seqTosample[s] = sampleName



derepInfoLevel_II = dict()


for line in open(sys.argv[2]):
    line = line.rstrip()
    data = line.split("\t")

    if len(data) > 1:
        derepInfoLevel_II[data[0].split("_")[0]] = [data[0]]+data[1].split()
    else:
        derepInfoLevel_II[data[0].split("_")[0]] = [data[0]]

outFile = open(sys.argv[4], 'w')

outFile.write("\t"+"\t".join(samples)+'\n')

for line in open(sys.argv[3]):
    data = line.split()
    counts = defaultdict(int)
    for d in data:
        for e in derepInfoLevel_II[d.split("_")[0]]:
            counts[seqTosample[e.split("_")[0]]] += int(e.split("_")[1])

    outFile.write("%s\t" % data[0].split("_")[0])
    for sample in samples:
        outFile.write("%s\t" % counts[sample])
    outFile.write("\n")

outFile.close()

