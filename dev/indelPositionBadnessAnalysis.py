import re
import matplotlib.pyplot as plt
import glob
import sys
# name match_code actual_table used_table orf match%

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


def main(input_dir, pat,):
    hits = "%s/%s" % (input_dir, pat)
    input_files = glob.glob(hits)
    max_len = 800
    graph_gap_score_by_position(input_files, '-', max_len)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: input_dir patt"
    else:
        main(*sys.argv[1:3])
        # bestHits = findBestHits("%s%s" % (data_dir, "/data/test2"))



                # take in a mhap out file
# pull out indels
# for each indel:
    # log position
    # if error score = -1 if not +1
    # avg all the scores for each position
    #

# for each indel, mark its position