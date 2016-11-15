import os
import re
import time
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo as matlist
from random import uniform, randrange, sample, seed

nucleotides = ['A', 'T', 'C', 'G'] * 2
seed(int(round(time.time() * 1000)))

def make_random_frame_shift(seq):
    f_shift_len = randrange(1, 3)
    if coinflip():
        return  insert_random_frame_shift(seq, f_shift_len)
    else:
        return delete_random_frame_shift(seq, f_shift_len)


def insert_random_frame_shift(victim, f_shift_len):
    noise = sample(nucleotides, f_shift_len)
    pos = rand_position(len(victim.seq))
    id = "%s_i%d@%d" % (victim.id, f_shift_len, pos)
    seq = str(victim.seq)[:pos] + "".join(noise) + str(victim.seq)[pos:]
    return SeqRecord(Seq(seq), id=id, description="")


def delete_random_frame_shift(victim, f_shift_len):
    pos = rand_position(len(victim.seq))
    seq = victim.seq[:pos] + victim.seq[(pos + f_shift_len):]
    id = "%s_d%d@%d" % (victim.id, f_shift_len, pos)
    return SeqRecord(seq, id=id, description="")


def insert_frame_shift_at(victim, f_shift_len, pos):
    noise = sample(nucleotides, f_shift_len)
    return insert_frame_shift_k_at(victim, noise, pos)


def insert_frame_shift_k_at(victim, noise, pos):
    seq = str(victim.seq)
    bad_seq = seq[:pos] + "".join(noise) + seq[pos:]
    victim.seq = Seq(bad_seq)
    #id = "%s_i%d@%d" % (victim.id, len(noise), pos)
    #victim.id = id
    return victim

def delete_frame_shift_at(victim, f_shift_len, pos):
    seq = str(victim.seq)
    bad_seq = seq[:pos] + seq[(pos + f_shift_len):]
    victim.seq = Seq(bad_seq)
    #id = "%s_d%d@%d" % (victim.id, f_shift_len, pos)
    #victim.id = id
    return victim

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


def getGC_from_name(aa_name):
    return aa_name.split('_')[-1]


def globalNucAlign(seq1,seq2):
    return pairwise2.align.globalxx(seq1, seq2)


def globalProtAlign(seq1, seq2, gap_open_penalty=-2, gap_extend_penalty=-1):
    return pairwise2.align.globalds(seq1, seq2, matlist.blosum62, gap_open_penalty, gap_extend_penalty)[0]


def isIndel(query_name):
    return '@' in query_name and ('_i' in query_name or '_d' in query_name)


def getIndelPosition(query_name):
    reg= r'.*_[id]\d@(\d+)'
    position = re.search(reg , query_name)
    return int(position.group(1))


def compute_gap_score(ref_align, query_align, gap_char):
    ref_gaps = ref_align.rstrip(gap_char).count(gap_char)
    query_align.rstrip(gap_char)
    query_gaps = query_align.rstrip(gap_char).count(gap_char)
    return ref_gaps + query_gaps

def getFileName(path):
    """Returns the filename (no extension, no directory) in an absoloute filepath.

    :param path:The file path.
    :return:    Returns the filename on a path with no directory prefix or file extension.
    """
    return os.path.splitext(os.path.basename(path))[0]


def translate_to_prot_by_name(fna_record, names_dict):
    gc = getGC_from_name(names_dict[fna_record.id])
    prot_seq = fna_record.seq.translate(table=gc)
    fna_record.seq = prot_seq
    return fna_record


def parseNames(names_file):
    names = {}
    for line in open(names_file):
        name, alias = line.split()
        names[name.rstrip()] = alias.rstrip()
    return names