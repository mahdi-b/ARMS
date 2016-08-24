import sys
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from random import uniform, randrange, sample, seed

nucleotides = ['A', 'T', 'C', 'G'] * 2

def frame_shift(seq):
    f_shift_len = randrange(1, 3)
    if coinflip():
        return  insert_frame_shift(seq, f_shift_len)
    else:
        return delete_frame_shift(seq, f_shift_len)


def insert_frame_shift(victim, f_shift_len):
    noise = sample(nucleotides, f_shift_len)
    pos = rand_position(len(victim))
    id = "%s_i%d@%d" % (victim.id, f_shift_len, pos)
    seq = str(victim.seq)[:pos] + "".join(noise) + str(victim.seq)[pos:]
    return SeqRecord(Seq(seq), id=id)


def delete_frame_shift(victim, f_shift_len):
    pos = rand_position(len(victim))
    seq = victim.seq[:pos] + victim.seq[(pos + f_shift_len):]
    id = "%s_d%d@%d" % (victim.id, f_shift_len, pos)
    return SeqRecord(seq, id=id)


def make_chimera(victims, victim):
    pos = rand_position(len(victim.seq))
    sibling = victims[randrange(0, len(victims))]
    seq = victim.seq[:pos] + sibling.seq[pos:]
    id = "%s_c_%s@%d" % (victim.id, sibling.id, pos)
    return SeqRecord(seq, id=id)


def coinflip():
    return int(uniform(0,1) * 100) % 2


def rand_position(len):
    return randrange(0, len)


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
    victims = get_victims(input_fasta, num_seqs)
    bads = []
    count = 0
    seed(int(round(time.time() * 1000)))

    for victim in victims:
        if coinflip():
            bads.append(make_chimera(victims, victim))
            count += 1
        elif coinflip():
            bads.append(frame_shift(victim))
            count += 1
        else:
            bads.append(victim)
    SeqIO.write(bads, open(output_file, 'w'), "fasta")
    print "Made %d bad seqs." % count

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: input_fasta  num_seqs  output_file"
        exit()
    print sys.argv[2]
    make_bads(*sys.argv[1:3])