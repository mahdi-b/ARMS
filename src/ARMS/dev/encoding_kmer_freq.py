import operator
import sys
from Bio import SeqIO
from collections import defaultdict
from itertools import product

groups = ["HRK", "DENBZKQX", "C", "STPAG", "MIJLV", "FYW"]
letters = ["A",   "B",     "C",  "D",     "E",    "F"]
prot_vals = {"A":1, "B":3, "C":5, "D":7, "E":14, "F":16}


def make_exchange_group():
    translator = {}
    for i in range(len(groups)):
        for protien in groups[i]:
            translator[protien] = letters[i]
    keys = product(letters, letters)
    return translator

def uset(dictionary, k):
    gset = set()
    keys = dictionary.keys()
    keys.remove(k)
    kmers = []
    for key in keys:
        gset |= set([k for (k,v) in dictionary[key]])
    return gset


def encode_fasta(input_file, k):
    dictionary = make_exchange_group()
    data = {}
    # counts_keys = ["".join(x) for x in product(letters, repeat=k)]
    # print counts_keys
    counts = {}
    names = []
    j = 0
    for seq in SeqIO.parse(open(input_file, 'r'), 'fasta'):
        j+=1
        if j%1000==0:
            print "%d\n"%j

        names.append(seq.description)

        # exchange group

        #translation = "".join([dictionary[x] for x in str(seq.seq).replace('X','')])
        translation = seq.seq
        #print translation
        # translation = clean_seq(seq)

        #counts = dict.fromkeys(counts_keys, 0)
        counts = defaultdict(int)
        # do freq analysis of size 2
        pos = range(len(translation) - k)
        for i in pos:
            counts["".join(translation[i:i+k])] += 1
        rslt = [(kmer,freq) for (kmer,freq) in counts.items() if freq > 0]
        # =====================encode top freqs===============================
        data[seq.description] = [kmer for (kmer, freq) in rslt]
    #print data
    #
    phyla = defaultdict(list)
    classes = defaultdict(list)
    for name in data.keys():
        pname = name.split("\t")[1].split(",")[0].split(":")[1]
        cname = pname

        try:
            cname = name.split("\t")[1].split(",")[1].split(":")[1]
        except:
            pass
        phyla[pname].append(data[name])
        classes[cname].append(data[name])

    print phyla.keys()
    print classes.keys()

    tax_level = classes
    common = {}
    for tax in tax_level.keys():
        membership = len(tax_level[tax])
        if membership > 5:
            print "%s (%d)" % (tax, membership)
            total = defaultdict(int)
            for member in tax_level[tax]:
                for kmer in member:
                    total[kmer] += 1

            rslt = sorted([(key,val) for (key,val)  in total.iteritems() if val >= membership*.7], key=operator.itemgetter(0))
            print rslt
            common[tax] = rslt
            # print len(rslt)

    """
    # look for uniqe kmers
    for key in common.keys():
        print key
        print "%d - %d" % (len(common[key]) ,len(uset(common,key)))
        mykmers = set([k for (k,v) in common[key]]) - uset(common,key)
        if len(mykmers) ==0:
            print "nope."
            exit()
    """
    def getkmerset(d):
        return set([k for (k,v) in d])

    for key1 in common.keys():
        print common[key1]
        for key2 in common.keys():
            print common[key2]
            if key1 != key2:
                print "%s vs %s" % (key1, key2)
                print getkmerset(common[key1])
                print getkmerset(common[key2])
                diff = len(getkmerset(common[key1]) - getkmerset(common[key2]))
                print diff
                if diff == 0:
                    print "nope."
                    exit()
    """
    dat = DataFrame(data).transpose()
    dat = dat
    phylum = []
    classs = []
    for name in dat.index.values:
        name = name.split("\t")[1]
        phylum.append(name.split(",")[0].split(":")[1])
        try:
            classs.append(name.split(",")[1].split(":")[1])
        except:
            classs.append(name.split(",")[0].split(":")[1])

    dat.loc[:, 'phylum'] = Series(phylum, index=dat.index)
    dat.loc[:, 'classs'] = Series(classs, index=dat.index)
    # Mollusca   Arthropoda
    # reduced_dat = reduced_dat.loc[reduced_dat['phylum'] == "Mollusca"]
    # reduced_dat = reduced_dat.loc[reduced_dat['classs'] == "Gastropoda"]
    # print reduced_dat.sort_values(cols)
    phylum_set = set(dat["phylum"])
    classs_set = set(dat["classs"])

    for phylum in phylum_set:
        total = set()
        ptable = dat.loc[dat['phylum'] == phylum]
        #print ptable
        for row in ptable.iterrows():
            print row['1']
    """
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: input_file k"
    else:
        encode_fasta(sys.argv[1], k=int(sys.argv[2]))
