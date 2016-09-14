import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from collections import defaultdict
from operator import itemgetter
from itertools import product
from computeDist import *

# Variables should be local not global..
# indel penalty.
gap_open = -10
# mismatch penalty
gap_extend = -0.5

# Index reference and query fastas
data_dir = ""
refs = ""
queries = ""

matrix = matlist.blosum62
# We assume that the "*" in a the query is an error of tanslation
# that would have been conserved amino acid otherwise 
allAAs = set([x[0] for x in matrix.keys()])
matrix.update(dict([((x, "*"), matrix[(x,x)]) for x in allAAs]))
matrix.update(dict([(("*", "*"), matrix[(x,x)]) for x in allAAs]))


def make_freq_dict(mhap_hits, ref, ref_index, nuc_prt_names_map):
    """Given a list of mhap_output hits for a query sequence, constructs a frequency table (dict) for the GC of those
        hits.
    :param mhap_hits:   A list of hits (which are tuples of (id, simmilarity) )
    :param ref:    A SeqIO.index of the reference nucleotide DB.
    :param ref_index:    A list of the keys in ref
    :return: An ordered list of tables to try
    """
    table_counts = {}
    for hit in mhap_hits:
        match_id = hit[0]
        match_nuc_name = returnRefById(match_id, ref_index, ref).id
        match_prot_name = nuc_prt_names_map[match_nuc_name]
        match_table = match_prot_name.split("_cd_")[1]

        if match_table in table_counts:
            table_counts[match_table] += 1
        else:
            table_counts[match_table] = 1
    return table_counts

def guess_table_by_freq(mhap_out_hits, ref, ref_index, nuc_prt_names_map, query_seq):
    """Given a list of mhap_output hits for a query sequence, and the query SeqRecord, guess the table number for that
        sequence based on the most commonly occuring GC in its mhap hits.

    :param mhap_out_hits:   A dictionary mapping sequence ID to a list of hits (which are tuples of (id, simmilarity) )
    :param ref:    A SeqIO.index of the reference nucleotide DB.
    :param ref_index:    A list of the keys in ref
    :param query_seq:   A query SeqRecord
    :return: An ordered list of tables to try
    """
    hits = mhap_out_hits[query_seq.id]
    table_counts = make_freq_dict(hits, ref, ref_index, nuc_prt_names_map)
    rslts = sorted(table_counts, key=table_counts.get, reverse=True)
    return ([int(x) for x in rslts], table_counts)


def guess_table_by_best_hit(mhap_out_hits, ref, ref_index, nuc_prt_names_map, query_seq):
    """Given a lsit of mhap_output hits for a query sequence, and a sequence, guess the table number for that sequence.

    :param mhap_out_hits:   A dictionary mapping sequence ID to a list of hits (which are tuples of (id, simmilarity) )
    :param ref:    A SeqIO.index of the reference nucleotide DB.
    :param ref_index:    A list of the keys in ref
    :param query_seq:   A query SeqRecord
    :return: An ordered list of tables to try
    """
    hits = mhap_out_hits[query_seq.id]
    best_hit = hits[0]
    for hit in hits:
        if hit[1] > best_hit[1]:
            best_hit = hit

    match_id = best_hit[0]
    match_nuc_name = returnRefById(match_id, ref_index, ref).id
    match_prot_name = nuc_prt_names_map[match_nuc_name]
    match_table = match_prot_name.split("_cd_")[1]

    table_counts = make_freq_dict(hits, ref, ref_index, nuc_prt_names_map)

    return ([int(match_table)], table_counts)





def mapAllHits(hitsFile):
    hits = {}
    for line in open(hitsFile):
        data = line.split()
        id = data [0]
        item = (data[1], float(data[3]))
        # if match in dict: append
        if id in hits.keys():
            hits[id].append(item)
        # if not in dict: install
        else:
            hits[id] = [item]

    return hits


def findBestHits(hitsFile):
    """Finds the best match score for each sequence id against all matches in the MHAP out filee

    :param hitsFile: A MHAP out put file
    :return: A dictionary containing the ID of the best match for each sequence in the MHAP result file.
    """
    bestHits=defaultdict(lambda: [0,0])
    # make sure file exists
    for line in open(hitsFile):
        data = line.split()
        # if current data[3]
        if bestHits[data[0]][1] < float(data[3]):
            bestHits[data[0]] = ( data[1], float(data[3]))
    return bestHits


def returnRefById(seqId, refsPositionsInFile, refs):
    """Returns a SeqRecord corresponding to the given sequence Id (ex. 1,2,3,.... n) from the reference fasta file.
    :param seqId: Sequence name as a string
    :param refsPositionsInFile: MAHP database ID
    :param refs: The SeqIO.index object of the reference database
    :return: The SeqRecord corresponding to that entry
    """
    refId  = refsPositionsInFile[int(seqId)-1]
    return refs[refId]

def returnQueryById(queryId, queries, orientation=0):
    """Returns the correctly oriented SeqRecord corresponding to the sequence named in the query fasta.
    :param queryId: Sequence name
    :param queries: The SeqIO.index object of the query fasta
    :param orientation: 1 or 0 representing orientation, 0 for forward, 1 for reverse
    :return: The correctly oriented SeqRecord from the query fasta.
    """
    # queryId is the true reference sequence
    if orientation == 1:
        return  queries[queryId].reverse_complement()
    else:
        return queries[queryId]


def globallyAlign(seq1, seq2, matrix=matrix, gap_open=gap_open, gap_extend=gap_extend):
    """Performs a global alignment for two sequences, and returns a score for simmilarity.

    :param seq1: String representation of a sequence
    :param seq2: String representation of a sequence
    :param matrix: Amino acid simmilarity scoring matrix
    :param gap_open: Integer score penalty for an open gap
    :param gap_extend: Integer score penalty for mismatch
    :return:
    """
    # seq1 and Seq2 are both strings, not Seq or SeqRecrod objects
    # get the aligned sequence
    ali = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)[0]
    # TODO add sim score here as was computed by Greg
    sim = (sum([1 for x in range(len(ali[0])) if ali[0][x] == ali[1][x] ]) * 1.0) / len(ali[0])
    return (ali, sim)



def run(nuc_ref_file, query_file, mhap_out_file, nuc_to_prot_file, prot_ref_file, heuristic_option):
    """

    :param nuc_ref_file: Nucleotide reference fasta
    :param query_file:  Query fasta
    :param mhap_out_file: MHAP output from running over the query_object file
    :param nuc_to_prot_file: A tab delimited mapping between nucleotide names in nuc_ref_file and protien names in prot_ref_fasta
    :param prot_refs: Amino acid refrence fasta
    :param heuristic_option: Int indicating which heuristic (if any) to use.  0: no heuristic, use tables 5,9,6
                                    1: use most frequently occuring GC in mhap hits, 2: use the GC of the best mhap hit.
    :return:
    """
    heuristic_option = int(heuristic_option)
    print "Parsing data...."
    nuc_refs = SeqIO.index(nuc_ref_file, 'fasta')
    nuc_refs_index = list(nuc_refs.keys())
    prot_refs = SeqIO.index(prot_ref_file, 'fasta')
    prot_refs_index = list(prot_refs.keys())
    queries = SeqIO.index(query_file, 'fasta')
    bestHits = findBestHits(mhap_out_file)
    print "Done parsing."

    # for heuristics
    h_name = "d"
    if heuristic_option > 0:
        h_name = "h"
    out = open("%s_%s.tab"% (h_name, mhap_out_file), 'w')
    alignment = ""
    allHits = mapAllHits(mhap_out_file)
    nuc_to_prot_map = parseNames(nuc_to_prot_file)
    i = 0
    # predict the best ORF for each query_object sequence
    for queryId in bestHits:
        print "=" * 200
        i += 1
        # Get correctly oriented SeqRecord
        query_object = returnQueryById(queryId, queries)
        # Get reference sequence
        print "%d/%d : %s" % (i, len(bestHits), queryId)
        nuc_query_name = queryId.split('_')[0]
        correct_table = int(nuc_to_prot_map[nuc_query_name].split("_cd_")[1])
        used_table = -1
        all_hits_table_freqs = {}
        candidate_tables = [5, 9, 6]
        if heuristic_option == 1:
            candidate_tables, all_hits_table_freqs = guess_table_by_freq(allHits, nuc_refs, nuc_refs_index, nuc_to_prot_map,
                                                                         query_object)
        if heuristic_option == 2:
            candidate_tables, all_hits_table_freqs = guess_table_by_best_hit(allHits, nuc_refs, nuc_refs_index,
                                                                         nuc_to_prot_map,
                                                                         query_object)
        print "Candidate tables: %s" % str(candidate_tables)
        found_good_translation = False
        query_seq_translation = ""
        # Try a permutation of open reading frames and translation candidate_tables
        found_good_translation = False
        for orf, table in product(range(3), candidate_tables):
            query_seq_translation = str(query_object[orf:].seq.translate(table=table))
            # If no stop codons, we're done
            if query_seq_translation.count("*") == 0:
                # TODO add another condition that the translation needs to match alignment landmarks
                # by aligning the reads on the protein sequence
                # doing a global alignment here
                found_good_translation = True
                print "Used orf: %d, table %d" % (orf, table)
                used_table = table
                break
        # Print successful translation
        simmilarity = 0
        if found_good_translation:
            status_code = 1
            if int(correct_table) != int(used_table):
                print "WRONG_TABLE_DETECTED!! %s : used %d instead of %d" % (queryId, used_table, correct_table)
                status_code = -1

            best_hit_id = bestHits[queryId][0]
            best_hit_protien_obj = returnRefById(best_hit_id, prot_refs_index, prot_refs)
            print "%s vs. %s" % (query_object.id, best_hit_protien_obj.id)

            alignResults = globallyAlign(best_hit_protien_obj.seq, query_seq_translation)
            simmilarity = 1 - float(computeDist_outterGapConserved([str(alignResults[0][0]), str(alignResults[0][1])]))
            print "Simmilarity: %f" % simmilarity
            print "Alignment: "
            print alignResults[0][0]
            print alignResults[0][1]
            alignment = "%s\t%s" % (alignResults[0][0], alignResults[0][1])
        # Welp, we tried.
        else:
            status_code = 0
            #print the ref table
            print "Error no translation found for sequence %s" % query_object.id
            print "Actual table was: %d" % correct_table
            # name match_code actual_table used_table orf match%
        #out_tab_line = "%s\t%d\t%d\t%d\t%d\t%f\n" % (queryId, status_code, correct_table, used_table, orf, simmilarity)
        topHits = allHits[queryId]
        topHits.sort(key=lambda x: x[1], reverse=True)
        named_hits =[]
        for hit in topHits[:5]:
            num_matched_kmers = hit[1]
            nuc_name = returnRefById(hit[0], nuc_refs_index, nuc_refs).id
            prot_name = nuc_to_prot_map[nuc_name]
            named_hits.append((prot_name, num_matched_kmers))
        count_dict = sorted(all_hits_table_freqs.items(), key=itemgetter(1), reverse=True)
        out_tab_line = "%s\t%d\t%s\t%s\t%d\t%s\n" % (queryId, correct_table, str(count_dict), str(named_hits), status_code, alignment)
        out.write(out_tab_line)


if __name__ == "__main__":
    if len(sys.argv) < 7:
        print "Usage: nuc_ref_fasta  query_fasta  mhap_out_file nuc_to_prot_file  prot_ref_fasta  heuristic_T_F"
    else:
        run(*sys.argv[1:7])
    #bestHits = findBestHits("%s%s" % (data_dir, "/data/test2"))



