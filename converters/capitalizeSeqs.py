import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
def capitalizeSeqs(input_fasta, output_fasta):
    seqBuffer = []
    i = 0
    if os.path.isfile(output_fasta):
        os.remove(output_fasta)

    output = open(output_fasta,'a')

    for sequence in SeqIO.parse(open(input_fasta,'rU'), "fasta"):
        sequence.seq = Seq(str(sequence.seq).upper(), SingleLetterAlphabet())
        seqBuffer.append(sequence)
        i += 1
        if i % 5000 == 0:
            SeqIO.write(seqBuffer, output, "fasta")
            seqBuffer = []
    SeqIO.write(seqBuffer, output, "fasta")
    output.close()