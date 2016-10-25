from Bio import SeqIO
import sys

def uniq(input_file, output_file, format):
    seen = {}
    output = []
    i = 0
    with open(output_file, 'w') as out_file:
        for seq in SeqIO.parse(open(input_file, 'r'), format):
            i += 1
            if not seen.has_key(seq.id):
                seen[seq.id] = 1
                output.append(seq)
                if i % 10000 == 0:
                    SeqIO.write(output, out_file, format)
                    output = []
        SeqIO.write(output, out_file, format)

if len(sys.argv) < 4:
    print "Usage: input_file, output_file, format\n"
    exit()
uniq(sys.argv[1], sys.argv[2], sys.argv[3])