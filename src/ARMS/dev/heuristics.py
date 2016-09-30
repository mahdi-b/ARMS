import matplotlib.pyplot as plt
import operator
import sys
from Bio import SeqIO
from functools import partial
from itertools import product
from pylab import rcParams

from utils import *


def freq_analysis(input_faa, window_size, graph=False):
    freq = {}
    freq ["*"] = 0
    total_kmers_observed = 0
    for record in SeqIO.parse(open(input_faa, 'r'), "fasta"):
        ref = str(record.seq)
        for i in range(len(ref) - window_size):
            total_kmers_observed += 1
            tup = ref[i:i + window_size]
            if "*" in tup:
                freq["*"] += 1
            if tup in freq.keys():
                freq[tup] += 1
            else:
                freq[tup] = 1
    for key in freq.keys():
        count = float(freq[key])
        freq[key] = count/total_kmers_observed
    sorted_x = sorted(freq.items(), key=operator.itemgetter(1), reverse=True)
    if not graph:
        return sorted_x , total_kmers_observed
    rng = range(len(sorted_x))
    keys = [sorted_x[i][0] for i in rng]
    vals = [sorted_x[i][1] for i in rng]
    print sorted_x
    #print freq['TAA']
    #print freq['TAG']
    #print freq['TGA']
    print len(sorted_x)
    rcParams['figure.figsize'] = 10, 50
    plt.barh(rng, vals)
    plt.yticks(rng, keys, fontsize=1)
    plt.xlabel('Freq')
    plt.title('Kmer frequency')
    plt.show()


def apply_to_fasta(input_fasta, output_fasta, fn_to_apply):
    out = open(output_fasta, 'w')
    records = []
    print "Generating %s..." % output_fasta
    for record in SeqIO.parse(open(input_fasta, 'r'), "fasta"):
        records.append(fn_to_apply(record))
    print "Writing %s..." % output_fasta
    SeqIO.write(records, out,'fasta')


def make_indel_fastas(input_fna):
    files = []
    out_template = getFileName(input_fna) + "_%s%d.%s"
    nuc_to_prot_names_file = "/home/greg/ARMS/data/bold110k_name_pairs.txt"
    names_dict = parseNames(nuc_to_prot_names_file)
    translate = partial(translate_to_prot_by_name, names_dict=names_dict)
    for i in [1,2]:
        # make insert shifts
        output_fna_name = out_template % ("i", i, "fna")
        # insert_frame_shift_k_at(victim, noise, pos):
        insert = partial(insert_frame_shift_k_at, noise='A'*i, pos=0)
        apply_to_fasta(input_fna, output_fna_name, insert)

        # translate
        output_faa_name = out_template % ("i", i, "faa")
        apply_to_fasta(output_fna_name, output_faa_name, translate)
        files.append(output_faa_name)

        # make delete shifts
        output_fna_name = out_template % ("d", i, "fna")
        # delete_frame_shift_at(victim, f_shift_len, pos):
        delete = partial(delete_frame_shift_at, f_shift_len=i, pos=0)
        apply_to_fasta(input_fna, output_fna_name, delete)

        # translate
        output_faa_name = out_template % ("d", i, "faa")
        apply_to_fasta(output_fna_name, output_faa_name, translate)
        files.append(output_faa_name)


    return files


def tabulate_indel_files(input_files):
    with open("permuted_rslts.txt", 'w') as out:
        for input_file, window_size in product(input_files, [2,3]):
            freqs, total = freq_analysis(input_file, window_size)
            out.write("%d\t%s\t%d\t%s\n" % (window_size, input_file, total, str(freqs)))

if __name__ == "__main__":
    """
    if len(sys.argv) < 3:
        print "Usage: faa_fasta window_size"
    else:
        freq_analysis(sys.argv[1], int(sys.argv[2]))
    """
    files = make_indel_fastas(sys.argv[1])
    files.append(sys.argv[1])
    tabulate_indel_files(files)