from __future__ import division

import random
import sys
from Bio import SeqIO

import numpy as np
from numpy.linalg import norm
from scipy.stats import entropy

DNA = {"A" :0b00, "C" :0b01, "G" :0b10, "T" :0b11, "N" :0b00}

def compute_certainty(votes):
    fwd_votes = sum(x > 0 for x in votes)
    rev_votes = sum(x < 0 for x in votes)
    if fwd_votes > rev_votes:
        return 1, float( fwd_votes /(len(votes ) *1.0))
    else:
        return -1, float( rev_votes /(len(votes ) *1.0))


def JSD(P, Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return -0.5 * (entropy(_P, _M) + entropy(_Q, _M))


def SumOfDiff(P ,Q):
    return np.sum(np.abs( P -Q))



def computeVal(myKmer):
    # TODO: make sure it won't overflow
    """
        Returns int value of a kmer by computing the binary value of its
        characters
    """
    binKmer = 0b00
    for c in myKmer:
        binKmer = binKmer <<  2
        try:
            binKmer = binKmer | DNA[c]
        except:
            pass
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
    counts = np.zeros( 4 **kmerSize)
    # TODO: update vlaue by shifting (<<) one char at a time
    for start in range(len(sequence )- kmerSize + 1):
        counts[computeVal(sequence[start: start +kmerSize])] += 1

    return counts


def guessOrientation(query_sequence, names, refs, k):

    nbSuccessFor = 0
    nbSuccessRev = 0
    nbTrials = 0
    for name in names:
        nbTrials +=1
        S_vec = returnKmersVec(query_sequence, k)
        R_vec = returnKmersVec(refs[name], k)
        forScore = SumOfDiff(S_vec, R_vec)
        S_rev_vec = returnKmersVec(query_sequence.reverse_complement(), k)
        revScore = SumOfDiff(S_rev_vec, R_vec)

        if forScore < revScore:
            nbSuccessFor += 1
        else:
            nbSuccessRev += 1

        if nbTrials > 15 and nbSuccessFor / nbTrials > 0.8:
            print nbSuccessFor, nbSuccessRev
            return "F"
        elif nbTrials > 15 and nbSuccessRev / nbTrials > 0.8:
            print nbSuccessFor, nbSuccessRev
            return "R"
    print nbTrials, nbSuccessFor, nbSuccessRev
    return "unresolved"


def runTest(size, boldSeqs, refs, k):
    names = np.random.choice(refs.keys(), 100)
    # refsChordata = SeqIO.to_dict(SeqIO.parse(open("/tmp/chordata.fasta", 'r'), "fasta"))
    # namesChordata = np.random.choice(refsChordata.keys() , 100)
    bioCodeSeqs = SeqIO.to_dict(
        SeqIO.parse(open("/tmp/BIOCODE_071216_MACSE_RENAMED_CORRECTED_FILTERD.fna", 'r'), "fasta"))

    seqNames = np.random.choice(bioCodeSeqs.keys(), size)

    for seqName in seqNames:
        if random.randint(0, 1):
            print "+"
            if guessOrientation(bioCodeSeqs[seqName], names, refs, k) != "F":
                print "failed in forward"
                print seqName
                print names
                # if guessOrientation(boldSeqs[seqName], namesChordata, refsChordata, k) != "F":
                #     print "failed in Chordata too"
                # else:
                #     print "Found in Chordata"
        else:
            print "-"
            if guessOrientation(bioCodeSeqs[seqName].reverse_complement(), names, refs, k) != "R":
                print "failed in reverse"
                print seqName
                print names
                # if guessOrientation(boldSeqs[seqName].reverse_complement(), namesChordata, refsChordata, k) != "R":
                #     print "failed in Chordata too"
                # else:
                #     print "Found in Chordata"
        print "-" * 80


if __name__ == "__main__":
    if len(sys.argv) < 3:
        # print "Usage: fast, k"
        boldSeqs = SeqIO.to_dict(SeqIO.parse(open("/home/gburgess/ARMS/data/bold10k.fna", 'r'), "fasta"))
        refs = SeqIO.to_dict(SeqIO.parse(open("/home/gburgess/ARMS/data/bold10k.fna", 'r'), "fasta"))
        runTest(5000, boldSeqs, refs, 7)
    else:
        run(sys.argv[1], int(sys.argv[2]))