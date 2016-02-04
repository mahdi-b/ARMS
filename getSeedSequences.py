from collections import defaultdict
import sys

seeds = defaultdict(list)


for line in open(sys.argv[1], 'r'):
    data = line.split()
    if data[0] == "H":
        seeds[data[-1].replace(":", "_")].append(data[-2].replace(":", "_"))


with open(sys.argv[2], 'w')  as outFile:
    for key in seeds:
        outFile.write(key+"\t")
        outFile.write(" ".join(seeds[key]))
        outFile.write("\n")



