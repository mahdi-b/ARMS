from Bio import SeqIO

seeds ={}
for line in open("8_parse_data/ALL_PARSED.out", 'r'):
    line = line.rstrip()
    seeds[line.split("\t")[0]] = line.split("\t")[1]

for mySeq in SeqIO.parse("8_parse_data/MACSEOUT_MERGED.fasta", 'fasta'):
    print mySeq.id
    if seeds.has_key(mySeq.id):
        for elem in seeds[mySeq.id].split(" "):
            print elem


