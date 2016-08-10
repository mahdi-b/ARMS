from collections import defaultdict
import sys

# TODO rename, these arent fasta files
def getSeedSequences(input_fasta, output_fasta):
    seeds = defaultdict(list)
    for line in open(input_fasta, 'r'):
        data = line.split()
        if data[0] == "H":
            seeds[data[-1].replace(":", "_")].append(data[-2].replace(":", "_"))


    with open(output_fasta, 'w')  as outFile:
        for key in seeds:
            outFile.write(key+"\t")
            outFile.write(" ".join(seeds[key]))
            outFile.write("\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: fastaA fastaB outfile"
        exit()
    getSeedSequences(sys.argv[1:3])