import operator
import sys
from Bio import SeqIO
from countToDict import parseCountFileToDict


def formatWithSwarmCounts(input_fasta, uc_parsed_out_file, output_fasta):
    # Covnerts a fasta and a dereplicated global count file to a sorted, dereplicated, fasta named by abundance
    seeds = []
    seedSizes = parseCountFileToDict(uc_parsed_out_file)

    print "Indexing reads"
    reads = SeqIO.index(input_fasta, "fasta")
    print "Done indexing reads"

    print "Renaming sequences"
    for name, count in sorted(seedSizes.items(), key=operator.itemgetter(1), reverse=True):
        print count
        s = reads[name]
        s.id = "%s_%s" % (name, count)
        s.description = ""
        seeds.append(s)
        print s
        # write in chunks
        if len(seeds) == 500000:
            SeqIO.write(seeds, open(output_fasta, 'a'), "fasta")
            seeds =[]
    # write the rest of the chunk
    SeqIO.write(seeds, open(output_fasta, 'a'), "fasta")
    print "Done renaming sequences"

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: fasta count_file outfile"
        exit()
    formatWithSwarmCounts(sys.argv[1:4])
