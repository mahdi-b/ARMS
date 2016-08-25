from alignReadsProt import *
from collections import defaultdict
from itertools import product

def guess_table(mhap_out_hits, refs, nuc_prt_names_map, query_seq):
    """Given a lsit of mhap_output hits for a query sequence, and a sequence, guess the table number for that sequence.

    :param mhap_out_hits:   A dictionary mapping sequence ID to a list of hits (which are split lines as a list)
    :param refs:    A SeqIO.index of the reference DB.
    :param query_seq:   A query SeqRecord
    :return: An ordered list of tables to try
    """
    table_counts = defaultdict(0)
    hits = mhap_out_hits[query_seq.id]
    for hit in hits:
        match_id = hit[1]
        match_nuc_name = returnRefById(match_id, refs)
        match_prot_name = nuc_prt_names_map[match_nuc_name]
        match_table = match_prot_name.split("_cd_")[1]
        table_counts [match_table] += 1
    print sorted(table_counts, key=table_counts.get)
    return sorted(table_counts, key=table_counts.get)

def parseNames(names_file, delim=','):
    names = {}
    for line in open(names_file):
        name, alias = line.split(delim)
        names[name.rstrip()] = alias.rstrip()
    return names
