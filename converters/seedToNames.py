from Bio import SeqIO


def seedToNames(seed_fasta, names_file):
    with open(names_file,'w') as output:
        for sequence in SeqIO.parse(open(seed_fasta,'rU'), "fasta"):
            output.write("%s\t%s\n" % (sequence.id, sequence.id))
    output.close()