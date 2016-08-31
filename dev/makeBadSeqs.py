import sys
import time
from utils import *


def get_victims(input_fasta, num_seqs):
    i = 0
    victims = []
    for seq in SeqIO.parse(input_fasta, 'fasta'):
        if i >= num_seqs:
            break
        seq.description = ""
        victims.append(seq)
        i += 1
    return victims


def make_bads(input_fasta, num_seqs, output_file):
    victims = get_victims(input_fasta, int(num_seqs))
    bads = []
    count = 0
    seed(int(round(time.time() * 1000)))

    for victim in victims:
        if coinflip:
            bads.append(make_random_frame_shift(victim))
            count += 1
        # elif coinflip():
        #   bads.append(make_chimera(victims, victim))
        #   count += 1
        else:
            bads.append(victim)
    SeqIO.write(bads, open(output_file, 'w'), "fasta")
    print "Made %d bad seqs." % count

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: input_fasta  num_seqs  output_file"
        exit()
    print sys.argv[2]
    make_bads(*sys.argv[1:4])