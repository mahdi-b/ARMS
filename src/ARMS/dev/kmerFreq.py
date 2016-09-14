from random import sample, seed
from Bio import SeqIO
import sys
from scipy.stats import entropy
from numpy.linalg import norm
import numpy as np
from datetime import datetime

DNA = {"A":0b00, "C":0b01, "G":0b10, "T":0b11, "N":0b00}

def compute_certainty(votes):
    fwd_votes = sum(x > 0 for x in votes)
    rev_votes = sum(x < 0 for x in votes)
    if fwd_votes > rev_votes:
        return 1, float(fwd_votes/(len(votes)*1.0))
    else:
        return -1, float(rev_votes/(len(votes)*1.0))


def JSD(P, Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return -0.5 * (entropy(_P, _M) + entropy(_Q, _M))


def which_way(fwd_score, rev_score):
    """Takes in a forward number and a reverse numebr and returns the higher of the two.
    """
    if fwd_score > rev_score:
        return 1
    elif fwd_score < rev_score:
        return -1
    else:
        return 0

def computeVal(myKmer):
    # TODO: make sure it won't overflow
    """
        Returns int value of a kmer by computing the binary value of its
        characters
    """
    binKmer = 0b00
    for c in myKmer:
        binKmer = binKmer <<  2
        binKmer = binKmer | DNA[c]
    return binKmer

def returnKmersVec(query, kmerSize, reverse=False):
    """
    returns an array of size 4^kmerSize, with counts for each kmer
    each kmer is index using its binary value.
    ex. ACGTA = 00 01 10 11 00 = 108
    ex. AAA   = 00 00 00 = 0
    """
    sequence = query.seq
    if reverse:
        sequence = str(query.seq)[::-1]
    counts = np.zeros(4**kmerSize)
    # TODO: update vlaue by shifting (<<) one char at a time
    for start in range(len(sequence)- kmerSize + 1):
        counts[computeVal(sequence[start:start+kmerSize])] += 1

    return counts


def get_orientation(querySeqRecord, refDict, k, inital_sample_size=10, max_sample_size=100):
    print "Generating samples.."
    available_sample_names = sample(refDict, max_sample_size)
    available_samples = [refDict[sample_name] for sample_name in available_sample_names]
    samples = available_samples[:inital_sample_size]
    certainty = 0
    orientation = 0
    query_freq_fwd = returnKmersVec(querySeqRecord, k)
    query_freq_rev = returnKmersVec(querySeqRecord, k, True)
    sample_freqs = [returnKmersVec(seq, k) for seq in samples]
    print "Done generating samples..."
    # inital check
    fwd_scores = [JSD(query_freq_fwd, sample_freq) for sample_freq in sample_freqs]
    rev_scores = [JSD(query_freq_rev, sample_freq) for sample_freq in sample_freqs]
    print "Fwdscores:"
    print fwd_scores
    print "Revscores:"
    print rev_scores
    current_sample_size = inital_sample_size

    votes = ([which_way(fwd_scores[i], rev_scores[i]) for i in range(current_sample_size)])
    orientation, certainty = compute_certainty(votes)
    print "computing values..."
    # while unsure...
    while certainty < 0.9 and current_sample_size < max_sample_size:
        print "UNSURE, DOING MORE TESTING %d %f" % (orientation, certainty)
        # add one more sample
        samples.append(available_samples[current_sample_size])
        current_sample_size += 1
        # do freq analysis on the new sample
        sample_freqs.append(returnKmersVec(samples[-1], k))
        fwd_scores.append(JSD(query_freq_fwd, sample_freqs[-1]))
        rev_scores.append(JSD(query_freq_rev, sample_freqs[-1]))
        votes.append(which_way(fwd_scores[-1], rev_scores[-1]))
        orientation, certainty = compute_certainty(votes)

    ending_sample_size = len(samples)
    if ending_sample_size >= 100 and certainty < 0.9:
        print "Tried %d samples, only %f sure." % (ending_sample_size, certainty)
    else:
        print "Tried %d samples, %f confident." % (ending_sample_size, certainty)
    return certainty, orientation


def run(referenceFasta, k):
    seed(datetime.now())
    print "Parsing file..."
    refs = SeqIO.to_dict(SeqIO.parse(open(referenceFasta, 'r'), "fasta"))
    print "Done parsing file"
    names = sample(refs, 100)
    queries = [refs[name] for name in names]

    fwds = queries[:50]
    revs = queries[50:]
    for query in revs:
        seq = query.seq
        query.seq = str(seq)[::-1]
        rslts = get_orientation(query, refs, k)
        print rslts


    for query in fwds:
        rslts = get_orientation(query, refs, k)
        print rslts

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: fast, k"
    else:
        run(sys.argv[1], int(sys.argv[2]))