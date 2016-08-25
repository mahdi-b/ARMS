from Bio import SeqIO, AlignIO
from itertools import chain, izip, combinations

import logging
import numpy as np
import re


def replaceInnerGap(line):
    """Replaces the internal gaps "-" of a string by ":".
        line:	"" A sequence to edit.

        return	"" The string with inner gaps replaced by ":".
    """
    frstA = re.search(r"^-*[^-]", line).end() - 1
    lastA = re.search(r"[^-]-*$", line).start()
    return line[:frstA] + (line[frstA:lastA].replace("-", ":")) \
           + line[lastA:]


def computeDist_outterGapConserved(aligns):
    """Computes the distance between the sequences in a .aln file,
        counting padding gaps as wildcard matches.

        aligns	 	[""] List of seqs to align.

        returns 	The mismatch % between aligned sequences.
    """
    miss = 0
    succ = 0
    for i in range(len(aligns)):
        aligns[i] = replaceInnerGap(aligns[i])

    for col in izip(*aligns):
        colSet = set(col)
        setSize = len(colSet)
        if ":" in colSet or (setSize > 2) or (setSize == 2 and not "-" \
                in colSet):
            miss = miss + 1
        else:
            succ = succ + 1
    return miss / (miss + float(succ))


def computeDist(aligns, distFn):
    """Computes the distance between sequences in a .aln file, using the
		distFn specified.

	 	aligns: 	[] List of seqs to align""
	 	distFn: 	"" Uppercase ID for the dist function to use
	 	
	 	return		The distance as a % between sequences in an 
					alignFile
	"""
    distFns = {"OUTTER_GAP_CONSERVED": computeDist_outterGapConserved}

    if distFns.has_key(distFn):
        return distFns[distFn](aligns)
    else:
        raise NameError("Error: Invalid grading Fn.")


def computePairwiseDistStats(alignedSeqs, distFn):
    """Does pairwise subsampling on a list of sequences.  Prints the
		distnace matrix and computes summary stats. 

	 	alignedSeqs:	[] The sequences to compare.

		return		(,) A tuple of (min,max,avg,std) distances
					summarizing the table.
	"""
    sampleSize = len(alignedSeqs)
    distMatrix = np.zeros((sampleSize, sampleSize))
    valList = []
    for i, j in combinations(range(sampleSize), 2):
        distMatrix[i, j] = computeDist([alignedSeqs[i], alignedSeqs[j]],
                                       distFn)
        valList.append(distMatrix[i, j])

    distMatrix
    myMin = min(valList)
    myMax = max(valList)
    myAvg = np.mean(valList)
    myStd = np.std(valList)
    sol = {"min": myMin, "max": myMax, "avg": myAvg, "std": myStd}

    logging.debug("Pairwise Data:\n")
    for key in sol.keys():
        logging.debug("%s: %f\n" % (key, sol[key]))
    return sol
