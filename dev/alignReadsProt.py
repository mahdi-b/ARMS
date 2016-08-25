import sys
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
import os
from collections import defaultdict
from itertools import product
from heuristics import guess_table, parseNames
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

def mapAllHits(hitsFile):
    hits = {}
    for line in open(hitsFile):
        data = line.split()
        id = data [0]
        # if match in dict: append
        if id in hits.keys():
            hits[id] = [data[1:]]
        # if not in dict: install
        else:
            hits[id].append(data[1:])
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


def returnRefById(seqId, refs):
    """Returns a SeqRecord corresponding to the given sequence Id (ex. 1,2,3,.... n) from the reference fasta file.

    :param seqId: Sequence name as a string
    :param refs: The SeqIO.index object of the reference database
    :return: The SeqRecord corresponding to that entry
    """
    refsPositionsInFile = list(refs.keys())
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


def run(ref_file, query_file, mhap_out_file, nuc_to_prot_file, prot_ref_fasta):
    refs = SeqIO.index(ref_file, 'fasta')
    queries = SeqIO.index(query_file, 'fasta')
    bestHits = findBestHits(mhap_out_file)

    # for heuristics
    allHits = mapAllHits(mhap_out_file)
    nuc_to_prot_map = parseNames(nuc_to_prot_file)
    prots = prot_ref_fasta
    i = 0
    # predict the best ORF for each query sequence
    for queryId in bestHits:
        i += 1
        # Get correctly oriented sequence
        query = returnQueryById(queryId, queries)
        # Get reference sequence
        ref = returnRefById(bestHits[queryId][0], refs)

        orfs = range(3)
        # tables = [5, 9, 6]
        tables = guess_table(allHits, refs, query, nuc_to_prot_map)
        foundTranslation = False
        translation = ""
        # Try a permutation of open reading frames and translation tables
        for orf, table in product(orfs, tables):
            foundTranslation = False
            translation = str(query[orf:].seq.translate(table=table))
            # If no stop codons, we're done
            if translation.count("*") == 0:
                # TODO add another condition that the translation needs to match alignment landmarks
                # by aligning the reads on the protein sequence
                # doing a global alignment here
                foundTranslation = True
                break
        # Print successful translation
        if foundTranslation:
            alignResults = globallyAlign(translation, ref.seq)

            print "%d/%d : %s vs. %s" % (i, len(bestHits), query.id, ref.id),
            print "-" * 200
            print alignResults[0][0]
            print alignResults[0][1]
            print alignResults[1]
            print "-" * 200
        # Welp, we tried.
        else:
            print("Error no translation found for sequence %s" % query.id)


seqTranslations = defaultdict(list)
if __name__ == "__main__":
    if len(sys.argv) < 6:
        print "Usage: ref_fasta  query_fasta  mhap_out_file nuc_to_prot_file  prot_ref_fasta"
    else:
        run(*sys.argv[1:6])
    #bestHits = findBestHits("%s%s" % (data_dir, "/data/test2"))



