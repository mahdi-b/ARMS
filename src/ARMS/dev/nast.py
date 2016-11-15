from collections import defaultdict
from classes.Helpers import init_pool, strip_ixes, printVerbose
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from Bio import SeqIO
from itertools import product
import subprocess
import select

def locate_deltas(ref_msa_template, nast_ref, nast_query, priority=0):
    """Finds extra characters in the MSA template and the regapped pairwise reference, logging them as deltas.

    :param ref_msa_template: Seq. The untouched MSA reference sequence.
    :param nast_ref: string.  The regapped pairwise reference sequence.
    :param nast_query: string.  The regapped pairwise query sequence.
    :param priority: int.  An integer that denotes which deltas belong together.  Should increase with each subsequent call.
    :return: tuple.  A tuple of (a dictionary mapping delta starting positions to the sequences inserted there,
                the total number of insertions in this query)
    """
    # NOTE: nast_ref is longer than ref_msa_template.  OR ELSE
    local_insertions = defaultdict(list)
    if len(ref_msa_template) == len(nast_ref):
        return local_insertions, 0

    if len(ref_msa_template) > len(nast_ref):
        print "ERROR: NAST SHORTER THAN REF MSA"
        exit()

    template_cursor = 0
    nast_cursor = 0
    template_insertions = 0
    while template_cursor < len(ref_msa_template):
        if ref_msa_template[template_cursor] != nast_ref[nast_cursor]:
            local_insertions[template_cursor].append(nast_query[nast_cursor])
            template_insertions += 1
        else:
            template_cursor += 1
        nast_cursor += 1
    if template_insertions + len(ref_msa_template) < len(nast_ref):
        print "ERROR: NOT ALL GAPS FOUND!"
        print "len MSA: %d" % len(ref_msa_template)
        print "TEMPLATE INSERTIONS: %d" %template_insertions
        print "len nast_ref: %d" % len(nast_ref)
        print "M: " + ref_msa_template.seq
        print "N: " + nast_ref
        print "Q: " + nast_query
        exit()
    # collapse the character insertion lists at each location (key) into a single string for insertion
    for key in local_insertions.keys():
        chars = "".join(local_insertions[key])
        local_insertions[key] = chars

    return local_insertions, template_insertions

def update_global_deltas(local_deltas, global_deltas):
    """Updates the list of global deltas with a set of local deltas (by appending the new deltas).

    :param local_deltas: dict{int:string}. A dictionary of local deltas, mapping insertion locations to insertion strings.
    :param global_deltas: dict{int:string}. A dictionary of global deltas, mapping insertion locations to insertion strings
    """
    for key in local_deltas.keys():
        global_deltas[key].add(local_deltas[key])


def mask_deltas(delta_dict, sequence):
    """Masks (by converting to lowercase) the deltas contributed by a sequence for easy identification later on.

    :param delta_dict: dict{int:string}. A dictionary of insertion strings indexed by position.  Contains only insertions contributed by
                        the sequence in question.
    :param sequence: string. The sequence to mask.
    :return: string. A masked copy of the sequence.
    """
    temp = list(sequence)
    insertions_so_far = 0
    keys = sorted(delta_dict.keys())
    for key in keys:
        num_inserts_here = len(delta_dict[key])
        pos = key + insertions_so_far
        temp[pos:pos + num_inserts_here]= list("".join(temp[pos:pos+num_inserts_here]).lower())
        insertions_so_far += num_inserts_here
    return "".join(temp)



def apply_deltas(delta_dict, sequence, gap_char='-', gap_letter=False):
    """Applies deltas to a sequence in order to make it fit in a new MSA.  Inserts gap_chars at the appropriate
    locations.  If a masked string (representing an insertion caused by this sequence) is detected, it is changed to
    uppercase.  Does not modify the passed in sequence.

    :param delta_dict: dict{int:string} A dictionary containing ALL changes that should be applied to this sequence.
    :param sequence: string. The sequence to use as a basis for modification.
    :param gap_char: char. The gap character to insert.
    :param gap_letter: bool. A custom flag used to generate the ruler.  If true, will always insert strings instead of gaps
                        into the sequence.
    :return: A copy of the sequence with the deltas applied.
    """
    temp = list(sequence)
    insertions_so_far = 0
    for key in sorted(delta_dict.keys()):
        pos = key + insertions_so_far
        max_char_len = max([len(s) for s in delta_dict[key]])
        gap_template = gap_char * max_char_len
        found = False
        possible_insertions = sorted(delta_dict[key], reverse=True)
        for chars in possible_insertions:
            length = len(chars)
            template = list(chars)
            template[-1:-1] = ([gap_char] * (max_char_len - length))
            # Option for making the ruler
            # print "Looking at %s, looking for %s @ %s" % ("".join(temp[pos:pos + length ]), chars.lower(), pos)
            if gap_letter:
                temp[pos:pos + max_char_len] = list(possible_insertions[0])
                found = True
            # If we find a delta this sequence contributed, just replace it (no insertion)
            elif "".join(temp[pos:pos + length ]) == chars.lower():
                temp[pos:pos + length ] = list(chars)
                temp[pos+length:pos+length] = [gap_char] * (max_char_len - length)
                found = True
                #print "replaced %s with %s" % (temp[pos+length:pos+length], template)
                break
            else: pass
        if not found:
            temp[pos:pos] = gap_template
        insertions_so_far += max_char_len
        # always jump a head, even if we didn't actually insert.
        # Because all sequences should line up to the same indexes at the end, they should all line up at any given point.
        # so we increment this to keep every cursor at the identical index at each iteration.

    return "".join(temp)


def get_realignment_regions(cumulative_insertions):
    """Given a string with deltas applied, find the gapping regions, and return them as a list of tuples.

    :param cumulative_insertions: dict{int:string}. A dictionary of all the deltas to apply.
    :return: [(),]. A list of tuples representing regions where insertions occured.
    """
    template = ' ' * 2 * max(cumulative_insertions.keys())
    gap_template = apply_deltas(cumulative_insertions, template, '@')
    gap_segments = []
    last_char = ''
    start = -1
    end = -1
    for i in range(len(gap_template)):
        if gap_template[i] == '@':
            if last_char == '@':
                end = i
            else:
                start = i
                end = i
        else:
            if last_char == '@':
                gap_segments.append((start, end))
        last_char = gap_template[i]
    return gap_segments


def muscle_realign(infile, outfile):
    """Calls muscle for realignment on a single region.  Returns a list of Seqs representing the given realigned region.

    :param infile: string. The input file path.
    :param outfile: string. The output file path.
    :return: [Seq,]. A list of realigned seqs.
    """
    select.select([open(infile)],[],[])
    cmd = "muscle -in %s -out %s" % (infile, outfile)
    subprocess.check_call(cmd, shell=True)
    select.select([open(outfile)],[],[])
    return [sequence.seq for sequence in SeqIO.parse(open(outfile, 'r'), 'fasta')]


def realign(seqs, start, end, outfilename):
    """Calls muscle with regions between start:end in each seq, and calls muscle on those regions.
    Returns a replacement dictionary for input sequences.
    :param seqs: [string,]  A list of all the sequences in an MSA.  These are sequences that have alrady had global
                            deltas applied.
    :param start: The start index of the realignment region.
    :param end: The end index of the realignment region.
    :param outfilename: The file to write the realignment sequences to.
    :return: dict{string:string}  A map of unique subsequences to their realignments.
    """
    old_vals = map(str.lower, list(set([seq[start:end].replace('-','') for seq in seqs])))

    try:
        old_vals.remove('')
    except:
        pass
    #print old_vals
    with open(outfilename, 'w') as output:
        out = ""
        for i in range(len(old_vals)):
            out += ">s%d\n%s\n" % (i, old_vals[i])
            i += 1
            if i % 5000 == 0:
                output.write(out)
                out = []
        output.write(out)
    rslt_name = "%s.mus.aln" % outfilename
    new_vals = muscle_realign(outfilename, rslt_name)
    #print new_vals
    val_map = dict((str(key), val) for key, val in zip(old_vals, new_vals))
    gaps = '-' * (end - start)
    val_map[''] = gaps
    #print val_map.keys()
    return val_map



def get_best_hits_from_vsearch(input_fna, ref_fna, outdir):
    """Calls vsearch with an input fasta, and returns a dictionary mapping each sequence to its best hit. (subject to
        the ID threshold (70%) in vsearch (See ProgramRunnerCommands.ALIGN_VSEARCH).

    :param input_fna: string.  Filepath to the input fna fasta file.
    :param ref_fna: string. Filepath to the reference fna fasta file.
    :param outdir: string. Filepath to the output directory for the hits file.
    :return: {string:string} A dictionary mapping input sequence names to the best hit in the reference DB.
    """
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

    vsearch_output = "%s/%s.out" % (outdir, strip_ixes(input_fna))

    # Choose the best hit
    return best_hits_from_vsearch(vsearch_output)


def make_faa_gc_lookup(name_map):
    """Reads a fna:faa names map file and returns a dict.  Also returns an fna-name:gc map.

    :param name_map: Filepath to a two-column file mapping fna sequence names to faa sequence names.
    :return: ( {faa names:fna names}, {faa-name:genetic code} )
    """
    lookup = {}
    fna_faa ={}
    for line in open(name_map,'r'):
        data = line.split("\t")
        lookup[data[0].rstrip()]  = int(data[1].rstrip().split('_')[-1])
        fna_faa[data[0].rstrip()] = data[1].rstrip()
    return lookup, fna_faa


def translate(record, is_fwd_read, gc):
    """Translates fna to faa.  Permutes over all ORF using the provided GC.

    :param record: SeqRecord.  The sequence to translate
    :param is_fwd_read: bool True if forward read.  Flase if RC.
    :param gc: int The genetic code table ot use.
    :return: A list of possible translations
    """
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

    :param input_fna: string. filepath to the input fna fasta file.
    :param best_hits: dict{name:name} map of fna sequence name to name of the closest match in the DB.
    :param gc_lookup_map: dict{name,int} map of fna sequence names to gc #.
    :return: dict{name:Seq} Map of faa name to translated Seq.
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



def augment(base_seq, augment):
    """Fills in mising data (gaps) in a base_seq with data from an augmenting sequence.  Does not modify inputs.

    :param base_seq: string. The base sequence to use as a template.  Non-gap caracters in this sequence will not be replaced.
    :param augment: string. The augmenting sequence to take data from.
    :return: [char,] A copy of base_seq where gaps are filled in with data from augment.
    """
    base = list(base_seq)
    for i in range(len(base)):
        if base[i] == '-':
            base[i] = augment[i]
    return base


def nast_regap(ref_msa_template, pairwise_ref, pairwise_query):
    """Applies the reference MSA's gapping to the pairwise query and ref sequences.

    :param ref_msa_template: string. The untouched msa sequence.
    :param pairwise_ref: string. the pairwise reference sequence.
    :param pairwise_query: string. the pairwise query sequence.
    :return: (a copy of the regapped query sequence, a copy of the regapped reference sequence)
    """
    """Taken from mothur's nast.cpp"""
    pairwiseLength = len(pairwise_query)
    candPair = list(pairwise_query)
    tempPair = list(pairwise_ref)
    tempAln = list(ref_msa_template)

    fullAlignIndex = 0
    pairwiseAlignIndex = 0

    newTemplateAlign = []
    candAln = []

    while (tempAln[fullAlignIndex] == '-'):
        candAln += '-'  ## add the initial '-' and '.' to the candidate and template
        newTemplateAlign += tempAln[fullAlignIndex]
        fullAlignIndex += 1

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
