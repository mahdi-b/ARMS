import matplotlib.pyplot as plt
from Bio import SeqIO
from multiprocessing import Pool

from makeBadSeqs import *
from utils import *


# name match_code actual_table used_table orf match%



def graph_gap_score_by_position(input_files, gap_char='-', max_len=800):
    total_indels = 0
    gaps_by_pos = [0]*max_len
    occurances_at_pos = [1]*max_len
    for file_ in input_files:
        with open(file_, 'r') as out:
            next(out)
            for line in out:
                data = line.split('\t')
                query_name = data[0]
                match_code = int(data[4])
                if isIndel(query_name):
                    total_indels += 1
                    pos = getIndelPosition(query_name)
                    print pos
                    if match_code == 1:
                        occurances_at_pos[pos] += 1
                        ref_align = data[5]
                        query_align = data[6]
                        gaps_by_pos[pos] += compute_gap_score(ref_align, query_align, gap_char)
    # average gaps per occurance + 1
    rslt = [ float(int(gaps_by_pos[i])/int(occurances_at_pos[i])) for i in range(max_len)]
    print rslt
    print "Graphing %d samples" % sum(occurances_at_pos)
    x = range(max_len)
    line, = plt.plot(x, rslt, '--', linewidth=2)
    line2, = plt.plot(x, occurances_at_pos, '-', linewidth=2)
    plt.show()



def translate_align_and_score(seqs):
    seq1, seq2, table, i = seqs
    prot1 = Seq(seq1).translate(table=table)
    prot2 = Seq(seq2).translate(table=table)
    rslt = globalProtAlign(prot1, prot2)
    return compute_gap_score(rslt[0][0], rslt[0][1], '-')


def fakeSample(input_fasta, max_len):
    total_indels = 0
    gaps_by_pos = [0]*max_len
    occurances_at_pos = [1]*max_len
    bad_seq = ""
    pool = Pool(4)
    work = []
    table = 5
    for record in SeqIO.parse(open(input_fasta, 'r'), "fasta"):
        ref = str(record.seq)
        for i in range(100):
            seq = ref
            if coinflip():
                bad_seq = insert_frame_shift_at(seq, 1, i)
            else:
                bad_seq = delete_frame_shift_at(seq, 1, i)
            x = (seq, bad_seq, table, i)
            work.append(x)
    rslt = pool.map_async(translate_align_and_score, work).get(9999)
    print rslt
    """
    print "Graphing %d samples" % sum(occurances_at_pos)
    x = range(max_len)
    line, = plt.plot(x, rslt, '--', linewidth=2)
    line2, = plt.plot(x, occurances_at_pos, '-', linewidth=2)
    plt.show()
    """

# grab 100 seqs
# for each seq,
# insert a char at position i:n
# align the new seq against itself
# grade teh alignment
def main(input_file):
    max_len = 303
    #graph_gap_score_by_position(input_files, '-', max_len)

    fakeSample(input_file, max_len)
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: input_file"
    else:
        main(*sys.argv[1:2])
        # bestHits = findBestHits("%s%s" % (data_dir, "/data/test2"))



                # take in a mhap out file
# pull out indels
# for each indel:
    # log position
    # if error score = -1 if not +1
    # avg all the scores for each position
    #

# for each indel, mark its position