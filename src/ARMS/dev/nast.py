from collections import defaultdict
import operator


#TODO: account for inserting into dictionary earlier than othe rpositions
# Todo account for inserting into the same space (append instead of overwrite)
def locate_insertions(ref_msa_template, nast_ref, nast_query, global_insertion_locations=defaultdict(list), priority=0):
    # NOTE: nast_ref is longer than ref_msa_template.  OR ELSE
    local_insertions = defaultdict(list)
    if len(ref_msa_template) == len(nast_ref):
        return global_insertion_locations, local_insertions

    if len(ref_msa_template) > len(nast_ref):
        print "ERROR: NAST SHORTER THAN REF MSA"
        exit()

    #print "template:\t%s\n" % ref_msa_template
    #print "nast_regapped:\t%s\n" % nast_ref
    template_cursor = 0
    nast_cursor = 0
    template_insertions = 0

    while template_cursor < len(ref_msa_template):
        a = ref_msa_template[template_cursor]
        b = nast_ref[nast_cursor]
        if a != b:
            #print "%d,%d: %s != %s" % (template_cursor, nast_cursor, a ,b)
            loc = template_cursor #+ template_insertions
            global_insertion_locations[loc].append((1, nast_query[nast_cursor], priority))
            local_insertions[loc].append((1, nast_query[nast_cursor], priority))
            template_insertions += 1
            nast_cursor += 1
        else:
            template_cursor += 1
            nast_cursor += 1

    if template_insertions + len(ref_msa_template) != len(nast_ref):
        print "ERROR: NOT ALL GAPS FOUND!"
        exit()

    return global_insertion_locations, local_insertions


def resolve_priority(insertions):
    # {164: [(1, 'V', 0)],
    rslt = []
    ins = sorted(insertions.items(), key=lambda x: x[0])
    for loc, todo in ins:
        for count, char, priority in todo:
            rslt.append((loc,count, char))
    return rslt



def delete_gaps(gap_locs, sequence):
    temp = list(sequence)
    gaps = sorted(gap_locs.items(), key=lambda x: x[0], reverse=True)
    for loc, (count, char) in gaps:
        temp.pop(loc)
    return "".join(temp)


def insert_gaps_ref(priority_list, sequence, gap_char='-', gap_letter=False):
    temp = list(sequence)
    insertions_so_far = 0
    for (loc, count, char) in priority_list:
        if gap_letter: gap_char = char
        temp.insert(loc + insertions_so_far, gap_char)
        insertions_so_far += count
    return "".join(temp)


def insert_gaps_query(priority_list, sequence, gap_char='-', gap_letter=False):
    temp = list(sequence)
    insertions_so_far = 0
    for (loc, count, char) in priority_list:

        if temp[loc + insertions_so_far] != char :
            if gap_letter: gap_char = char
            temp.insert(loc + insertions_so_far, gap_char)
        insertions_so_far += count
    return "".join(temp)
"""
def insert_gaps(gap_locs, sequence, gap_char='-', gap_letter=False):


    temp = list(sequence)
    insertions_so_far = 0
    gaps = sorted(gap_locs.items(), key=lambda x: x[2]*100 + x[1])
    print gaps
    for loc, inserts in gaps:
        "======================================"
        for (count, char) in inserts:
            if temp[loc+insertions_so_far] != char :
                if gap_letter: gap_char = char
                temp.insert(loc, gap_char)
            insertions_so_far += count
        "====================================="
    return "".join(temp)
"""

def to_cig(seq):
    last_char = seq[0]
    count = 1
    out = ""
    for char in seq[1:]:
        if char != last_char:
            out+="%d%s" % (count, last_char)
            last_char = char
            count = 1
        else:
            count += 1
    out+="%d%s" % (count, last_char)
    return out



def nast_regap(ref_msa_template, pairwise_ref, pairwise_query):
    """Taken from mothur's nast.cpp"""
    pairwiseLength = len(pairwise_query)
    fullAlignLength = len(ref_msa_template)
    candPair = list(pairwise_query)
    tempPair = list(pairwise_ref)
    tempAln = list(ref_msa_template)

    fullAlignIndex = 0
    pairwiseAlignIndex = 0

    newTemplateAlign = []
    candAln = []

    while (tempAln[fullAlignIndex] == '.' or tempAln[fullAlignIndex] == '-'):
        candAln += '-'  ## add the initial '-' and '.' to the candidate and template
        newTemplateAlign += tempAln[fullAlignIndex]
        fullAlignIndex += 1

    lastLoop = ""
    while (pairwiseAlignIndex < pairwiseLength):
        if tempPair[pairwiseAlignIndex].isalpha() and tempAln[fullAlignIndex].isalpha() and candPair[
            pairwiseAlignIndex].isalpha():
            #  the template and candidate pairwise and template aligned have characters
            #	need to add character onto the candidatSeq.aligned sequence

            candAln += candPair[pairwiseAlignIndex]
            newTemplateAlign += tempPair[pairwiseAlignIndex]

            pairwiseAlignIndex += 1
            fullAlignIndex += 1

        elif tempPair[pairwiseAlignIndex].isalpha() and not tempAln[fullAlignIndex].isalpha() and candPair[
            pairwiseAlignIndex].isalpha():
            #	the template pairwise and candidate pairwise are characters and the template aligned is a gap
            #	need to insert gaps into the candidateSeq.aligned sequence

            candAln += '-'
            newTemplateAlign += '-'
            fullAlignIndex += 1

        elif not tempPair[pairwiseAlignIndex].isalpha() and tempAln[fullAlignIndex].isalpha() and candPair[
            pairwiseAlignIndex].isalpha():
            #  the template pairwise is a gap and the template aligned and pairwise sequences have characters
            #	this is the alpha scenario.  add character to the candidateSeq.aligned sequence without progressing
            #	further through the tempAln sequence.

            candAln += candPair[pairwiseAlignIndex]
            newTemplateAlign += '-'
            pairwiseAlignIndex += 1

        elif tempPair[pairwiseAlignIndex].isalpha() and tempAln[fullAlignIndex].isalpha() and not candPair[
            pairwiseAlignIndex].isalpha():
            #  the template pairwise and full alignment are characters and the candidate sequence has a gap
            #	should not be a big deal, just add the gap position to the candidateSeq.aligned sequence

            candAln += candPair[pairwiseAlignIndex]
            newTemplateAlign += tempAln[fullAlignIndex]  #
            fullAlignIndex += 1
            pairwiseAlignIndex += 1

        elif not tempPair[pairwiseAlignIndex].isalpha() and not tempAln[fullAlignIndex].isalpha() and candPair[
            pairwiseAlignIndex].isalpha():
            #	the template pairwise and aligned are gaps while the candidate pairwise has a character
            #	this would be an insertion, go ahead and add the character->seems to be the opposite of the alpha scenario

            candAln += candPair[pairwiseAlignIndex]
            newTemplateAlign += tempAln[fullAlignIndex]
            pairwiseAlignIndex += 1
            fullAlignIndex += 1

        elif tempPair[pairwiseAlignIndex].isalpha() and not tempAln[fullAlignIndex].isalpha() and not candPair[
            pairwiseAlignIndex].isalpha():
            #	template pairwise has a character, but its full aligned sequence and candidate sequence have gaps
            #	this would happen like we need to add a gap.  basically the opposite of the alpha situation

            newTemplateAlign += tempAln[fullAlignIndex]
            candAln += "-"
            fullAlignIndex += 1

        elif not tempPair[pairwiseAlignIndex].isalpha() and tempAln[fullAlignIndex].isalpha() and not candPair[
            pairwiseAlignIndex].isalpha():
            #	template and candidate pairwise are gaps and the template aligned is not a gap this should not be possible
            #	would skip the gaps and not progress through full alignment sequence
            #	not tested yet

            print "We're into D " + fullAlignIndex + " " + pairwiseAlignIndex
            pairwiseAlignIndex += 1

        else:
            #	everything has a gap - not possible
            #	not tested yet

            print "We're into F " + fullAlignIndex + " " + pairwiseAlignIndex
            pairwiseAlignIndex += 1
            fullAlignIndex += 1

    trailing_dashes = ['-'] * (fullAlignLength - len(candAln))
    trailing_dashes = tempAln[fullAlignIndex:]
    candAln += trailing_dashes
    newTemplateAlign += trailing_dashes
    #newTemplateAlign += tempAln[i]

    candAln = "".join(candAln).upper()  # everything is upper case
    return candAln, "".join(newTemplateAlign).upper()