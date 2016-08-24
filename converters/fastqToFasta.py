import sys
from Bio import SeqIO

def translateFastqToFasta(inputFastQ, outputFasta):
    with open(inputFastQ, 'r') as input:
        with open(outputFasta,'w') as output:
            for record in SeqIO.parse(open(inputFastQ, 'r'), "fastq"):
                output.write(record.format("fasta"))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: input_fastq_file output_fasta_file"
        exit()
    translateFastqToFasta(sys.argv[1], sys.argv[2])