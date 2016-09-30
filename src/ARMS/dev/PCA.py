import matplotlib
import sys
from Bio import SeqIO
from itertools import product
from random import choice

from pandas import DataFrame, Series
from sklearn import decomposition

groups = ["HRK", "DENBZKQX", "C", "STPAG", "MIJLV", "FYW"]
letters = ["A",   "B",     "C",  "D",     "E",    "F"]
letters = [x for x in "HRKDENBZKQXCSTPAGMIJLVFYW"]
letters = ["A","T","C","G"]
nucs = ["A","T","C","G"]
nuc_vals ={"A":1, "T":100, "C":1000, "G":10000}
prot_vals = {"A":1, "B":3, "C":5, "D":7, "E":14, "F":16}

"""
R.................A or G
Y.................C or T
S.................G or C
W.................A or T
K.................G or T
M.................A or C
B.................C or G or T
D.................A or G or T
H.................A or C or T
V.................A or C or G
N.................any
"""

iupac_map = {"Y":["C","T"],
             "R":["A","G"],
             "S": ["C", "G"],
             "W": ["A", "T"],
             "K": ["T", "G"],
             "M": ["A", "C"],
             "B": ["C", "G", "T"],
             "D": ["A", "G", "T"],
             "V": ["C", "G", "A"],
             "H": ["A", "C", "T"],
             "N": ["A", "C", "T", "G"],
             }

def numeralize(seq, count, exchange=True):
    s=0
    for nuc in seq:
        #s *= 10
        if exchange:
            s += prot_vals[nuc]
        else:
            s += ord(nuc) - ord('A')
    return s


def make_exchange_group():
    translator = {}
    for i in range(len(groups)):
        for protien in groups[i]:
            translator[protien] = letters[i]
    keys = product(letters, letters)
    return translator


def clean_seq(seq):
    s = [x for x in seq]
    for i in range(len(s)):
        c = s[i]
        if c in nucs:
            pass
        else:
            s[i] = choice(iupac_map[c])
    return "".join(s)


def encode_fasta(input_file, k):
    #dictionary = make_exchange_group()
    data = {}
    counts_keys = ["".join(x) for x in product(letters, repeat=k)]
    # print counts_keys
    counts = {}
    names = []
    j = 0
    for seq in SeqIO.parse(open(input_file, 'r'), 'fasta'):
        j+=1
        if j%1000==0:
            print "%d\n"%j
        # translate to exchange group
        name = seq.description.split("\t")
        names.append(name[1])

        # exchange group

        #translation = "".join([dictionary[x] for x in str(seq.seq).replace('X','')])
        #translation = seq.seq
        translation = clean_seq(seq)

        counts = dict.fromkeys(counts_keys, 0)

        # do freq analysis of size 2
        pos = range(len(translation) - k)
        for i in pos:
            counts["".join(translation[i:i+k])] += 1
        # =====================encode top freqs===============================
        #data["".join(name[0].split(";")[:1])] = [numeralize(kmer,freq) for (kmer, freq) in sorted(counts.items(), key=operator.itemgetter(1), reverse=True)[:10]]


        # =========================encode full freq analysis===================
        data["".join(name[0].split(";")[:1])] = [counts[x] for x in counts_keys]
    return data, names


def do_pca(input_file, k, dim, graph=False):
    data, names = encode_fasta(input_file, k)
    dat = DataFrame(data).transpose()
    print dat.to_string()
    pca = decomposition.PCA()
    pca.fit(dat)
    vars = pca.explained_variance_
    vsum = sum(vars)
    normd = [i/vsum for i in vars]
    print normd
    print normd[:dim]
    print sum(normd[:dim])
    cols = ['x','y','z'][:dim]
    pca.n_components = dim
    reduced = pca.fit_transform(dat)
    reduced_dat = DataFrame(reduced, index=names, columns=cols)

    phylum =[]
    classs = []
    for name in reduced_dat.index.values:
        phylum.append(name.split(",")[0].split(":")[1])
        try:
            classs.append(name.split(",")[1].split(":")[1])
        except:
            classs.append(name.split(",")[0].split(":")[1])

    reduced_dat.loc[:, 'phylum'] = Series(phylum, index=reduced_dat.index)
    reduced_dat.loc[:, 'classs'] = Series(classs, index=reduced_dat.index)
    # Mollusca   Arthropoda
    #reduced_dat = reduced_dat.loc[reduced_dat['phylum'] == "Arthropoda"]
    #Gastropoda
    # reduced_dat = reduced_dat.loc[reduced_dat['classs'] == "Insecta"]
    # print reduced_dat.sort_values(cols)
    reduced_dat.to_csv(path_or_buf="%s.dataframe" % input_file.split(".")[0], sep="\t", columns= ["phylum", "classs"]+cols,
                       index=True, header=True)
    phylum_set = set(reduced_dat["phylum"])
    classs_set = set(reduced_dat["classs"])
    print classs_set
    n_tax = len(phylum_set)
    n_class = len(classs_set)
    colors =  [name for name, hex in matplotlib.colors.cnames.iteritems()]
    col_map = {key: value for (key, value) in zip(classs_set, colors[:n_class])}
    # Color the dataframe
    coloring = [col_map[classs] for classs in reduced_dat["classs"]]
    print phylum_set
    print n_tax
    if graph:
        #=================
        #== Plot 2d ======
        #=================
        if dim == 2:
            import seaborn as sns
            print "Graphing..."
            ax = sns.lmplot(x='x', y='y',data=reduced_dat, hue="classs", fit_reg=False)
            sns.plt.savefig("myfig.png")
            sns.plt.show()
        #=================
        #== Plot 3d ======
        #+================
        if dim == 3:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(reduced_dat['x'], reduced_dat['y'], reduced_dat['z'], c=coloring)
            plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: input_file k dim(2|3)"
    else:
        do_pca(sys.argv[1], k=int(sys.argv[2]), dim=int(sys.argv[3]),graph=True)

"""
1.use the groups mahdi gave, mapping protiens to A,B,C,D...
2. translate protien seq to the mapping
3. do kmer freq analysis of the seq
4 thats your vector
"""