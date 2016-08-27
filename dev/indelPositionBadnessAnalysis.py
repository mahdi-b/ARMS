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


def main(input_dir, pat):
    total_indels = 0
    max_len = 1000
    errors_by_pos = [0]*max_len
    totals_by_pos = [0]*max_len
    hits = "%s/%s" % (input_dir, pat)
    input_files = glob.glob(hits)
    for file_ in input_files:
        with open(file_, 'r') as out:
            next(out)
            for line in out:
                data = line.split()
                query_name = data[0]
                match_code = data[4]
                if isIndel(query_name):
                    print query_name
                    total_indels += 1
                    pos = getIndelPosition(query_name)
                    print pos
                    totals_by_pos[pos] += 1
                    if match_code in [0,-1]:
                        # found a stop
                        errors_by_pos[pos] += 1

    x = range(max_len)
    print "Bad Errors = %s/%s" % (sum(errors_by_pos) , total_indels)
    line, = plt.plot(x, errors_by_pos, '--', linewidth=2)
    plt.show()

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