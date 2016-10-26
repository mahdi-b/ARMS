from collections import defaultdict
from classes.Helpers import init_pool, strip_ixes, printVerbose
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from Bio import SeqIO
from itertools import product

def locate_deltas(ref_msa_template, nast_ref, nast_query, global_insertion_locations=defaultdict(list), priority=0):
    # NOTE: nast_ref is longer than ref_msa_template.  OR ELSE
    local_insertions = defaultdict(list)
    if len(ref_msa_template) == len(nast_ref):
        return global_insertion_locations, local_insertions

    if len(ref_msa_template) > len(nast_ref):
        print "ERROR: NAST SHORTER THAN REF MSA"
        exit()

    template_cursor = 0
    nast_cursor = 0
    template_insertions = 0
    order = 0
    while template_cursor < len(ref_msa_template):
        if ref_msa_template[template_cursor] != nast_ref[nast_cursor]:
            global_insertion_locations[template_cursor].append((template_cursor, nast_query[nast_cursor], priority+order))
            local_insertions[template_cursor].append((template_cursor, nast_query[nast_cursor], priority+order))
            template_insertions += 1
        else:
            template_cursor += 1
        nast_cursor += 1
        order += .001
    if template_insertions + len(ref_msa_template) != len(nast_ref):
        print "ERROR: NOT ALL GAPS FOUND!"
        exit()

    return global_insertion_locations, local_insertions


def resolve_deltas(insertions):
    for loc in insertions.keys():
        delta = sorted(insertions[loc], key=lambda x: x[2])
        insertions[loc] = delta
    return insertions


def mask_deltas(delta_dict, sequence):
    temp = list(sequence)
    insertions_so_far = 0
    for key in sorted(delta_dict.keys()):
        num_inserts_here = len(delta_dict[key])
        pos = key + insertions_so_far
        temp[pos:pos + num_inserts_here]= list("".join(temp[pos:pos+num_inserts_here]).lower())
        insertions_so_far += num_inserts_here
    return "".join(temp)


def apply_deltas(delta_dict, sequence, gap_char='-', gap_letter=False):

    temp = list(sequence)
    insertions_so_far = 0
    for key in sorted(delta_dict.keys()):
        rep_offset = 0
        pos = key + insertions_so_far
        for (loc, char, priority) in delta_dict[key]:
            # Option for making the ruler
            if gap_letter:
                temp[pos + rep_offset] = char
                rep_offset += 1
            # If we find a delta we contributed, just change it
            elif temp[pos + rep_offset] == char.lower():
                temp[pos+rep_offset] = char
                rep_offset += 1
            # Insert someone else's delta
            else: temp.insert(pos+ rep_offset, gap_char)
            # always jump a head, even if we didn't actually insert.
            # Because all sequences should line up to the same indexes at the end, they should all line up at any given point.
            # so we increment this to keep every cursor at the identical index at each iteration.
            insertions_so_far += 1

    return "".join(temp)


def get_best_hits_from_vsearch(input_fna, ref_fna, outdir):

    def best_hits_from_vsearch(v_search_output):
        best_hits = {}
        for line in open(v_search_output, 'r'):
            data = line.split("\t")
            query_name = data[0].rstrip()
            if best_hits.has_key(query_name):
                if float(best_hits[query_name][2].rstrip()) < float(data[2].rstrip()):
                    best_hits[query_name] = data
            else:
                best_hits[query_name] = data
        return best_hits


    threads = 1
    pool = init_pool(threads)
    #printVerbose.VERBOSE = True
    print "calling vsearch"
    processes=1
    aln_user_string=""
    extraargstring=""


    printVerbose("Aligning against reference sequences...")
    #     # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
    # --userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt
    ProgramRunner(ProgramRunnerCommands.ALIGN_VSEARCH,
                            [processes, input_fna, ref_fna, "%s/%s.out" % (outdir, strip_ixes(input_fna)),
                             "%s/%s.alnout" % (outdir, strip_ixes(input_fna)), aln_user_string],
                            {"exists": [input_fna, ref_fna], "positive": [processes]},
                            extraargstring).run()


    print "cleaning up."
    vsearch_output = "%s/%s.out" % (outdir, strip_ixes(input_fna))

    # Choose the best hit
    return best_hits_from_vsearch(vsearch_output)


def make_faa_gc_lookup(name_map):
    lookup = {}
    fna_faa ={}
    for line in open(name_map,'r'):
        data = line.split("\t")
        lookup[data[0].rstrip()]  = int(data[1].rstrip().split('_')[-1])
        fna_faa[data[0].rstrip()] = data[1].rstrip()
    return lookup, fna_faa


def translate(record, is_fwd_read, gc):
    possible_translations = []
    for orf, table in product(range(3), [gc]):
        translation = ""
        if is_fwd_read:
            translation = str(record[orf:].seq.translate(table=table))

        else:
            translation = str(record[orf:].reverse_complement().seq.translate(table=table))
        # If no stop codons, we're done
        if str(translation).count("*") == 0:
            possible_translations.append(translation)
    return possible_translations


def translate_to_prot(input_fna, best_hits, gc_lookup_map):
    """Returns dict of non-stopping translations for sequences in best_hits by name.

    :param input_fna:
    :param best_hits:
    :param name_map:
    :return:
    """
    translations = {}

    # grab all input seqs as a dict
    input_seqs = SeqIO.to_dict(SeqIO.parse(open(input_fna, "r"), "fasta"))

    for key in best_hits.keys():
        # BALI4606_0_ID5_656	BMOO-17606	98.1	311	100.0
        data = best_hits[key]
        query_name = data[0].rstrip()
        ref_hit_name = data[1].rstrip()
        is_fwd_read = data[-1].rstrip() == "+"
        gc = gc_lookup_map[ref_hit_name]
        myseq = input_seqs[query_name]
        possible_translations = translate(myseq, is_fwd_read, gc)
        translations[query_name] = possible_translations
    return translations


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

    trailing_dashes = tempAln[fullAlignIndex:]
    candAln += trailing_dashes
    newTemplateAlign += trailing_dashes

    return "".join(candAln).upper(), "".join(newTemplateAlign).upper()