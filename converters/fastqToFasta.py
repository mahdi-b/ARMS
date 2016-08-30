import sys
from Bio import SeqIO

def translateFastqToFasta(inputFastQ, outputFasta):
    """Converts a FASTQ file into a FASTA file.

    :param inputFastQ: The filepath to the input fastq file.
    :param outputFasta: The filepath to the output fasta file.
    :return: The filepath to the output fasta file.
    """
    records = []
    i = 0
    with open(inputFastQ, 'r') as input:
        with open(outputFasta,'w') as output:
            for record in SeqIO.parse(open(inputFastQ, 'r'), "fastq"):
                records.append(record)
                i += 1
                if i % 5000 == 0:
                    SeqIO.write(records, output, 'fasta')
                    records = []
            SeqIO.write(records, output, 'fasta')
    return outputFasta


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: input_fastq_file output_fasta_file"
    else:
        translateFastqToFasta(*sys.argv[1:3])