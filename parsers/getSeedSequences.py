from collections import defaultdict
import sys

# TODO rename, these arent fasta files
def getSeedSequences(input_uc_file, output_names_file):
    """Pulls the representative seed sequences from a uc file and generates a .names file.

    :param input_uc_file:
    :param output_names_file:
    :return:
    """

    seeds = defaultdict(list)
    for line in open(input_uc_file, 'r'):
        data = line.split()
        if data[0] == "H":
            seeds[data[-1].replace(":", "_")].append(data[-2].replace(":", "_"))


    with open(output_names_file, 'w')  as outFile:
        for key in seeds:
            outFile.write(key+"\t")
            outFile.write(" ".join(seeds[key]))
            outFile.write("\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: fastaA fastaB outfile"
    else:
        getSeedSequences(*sys.argv[1:3])