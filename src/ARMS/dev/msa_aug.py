from Bio import SeqIO
from collections import defaultdict
from consensus import dumb_consensus
from nast import nast_regap, locate_deltas, apply_deltas, mask_deltas, make_faa_gc_lookup, \
    get_best_hits_from_vsearch, translate_to_prot, update_global_deltas, get_realignment_regions, augment, realign
from utils import globalProtAlign

import os

#query.fna", bold10k.fna", bold10k.faa.msa, bold10k_name_pairs.txt
def add_to_msa(input_fna, ref_fna, ref_faa_msa, name_map, outdir):

    load_saved_data = True
    save_data_to_file = False

    # read MSA file to dict
    ref_msa = SeqIO.to_dict(SeqIO.parse(open(ref_faa_msa, 'r'), 'fasta'))

    # store the gc for each ref seq
    gc_lookup_map, fna_faa_map = make_faa_gc_lookup(name_map)

    # Find closest neighbors for each input sequence
    best_hits = get_best_hits_from_vsearch(input_fna, ref_fna, outdir)
    if len(best_hits) == 0:
        print "\n\n\nERROR: no matching IDs found in reference DB."
        exit()

    # Translate to Protiens
    translations = translate_to_prot(input_fna, best_hits, gc_lookup_map)
    name_keys = translations.keys()

    # make a list of msa ref seqs
    msa_template_refs = dict((name, ref_msa[fna_faa_map[best_hits[name][1]]]) for name in name_keys)

    # Load files if possible
    loaded = False
    savefilepath = "%s.save" % input_fna
    if os.path.isfile(savefilepath) and load_saved_data:
        lines = "".join([line.rstrip() for line in open(savefilepath, 'r')])
        loaded = True
        data = eval(lines)
        pairwise_queries = data["pairwise_queries"]
        pairwise_refs = data["pairwise_refs"]
        consensus = data["consensus"]
        t_consensus = list(consensus)
        aug_msa_template_refs = dict((name, augment(msa_template_refs[name], t_consensus)) for name in name_keys)

    else:
        print "Generating a consensus.."
        # Generate a consensus of the MSA
        # consensus = "-------------------------------------G---MMVMFADRWLFSTNHKDIGTLYFIFGAWIAGM--IVGTSLSL------------------LIRAELGQPPG-SHLLNI----GDNDQ--S--FL--------------IYNVIVTAHAFIMIFF--MVMPIMIGGFG--L-------NWLVPLMLNTSFFDPAGGAPDMAFPRMNNMS---F---K--FW--LLGHPPSLTLLLSSS-----M-----I-VENGRKK-F--KC-----TA-------KKEAFGTGWTVYP---PLSSNILGFVVWAHSGASVD----------------DLAI------FS----------------------LHLAGISSILG---------------------------AINFI--------------TTIINMSR--IISNGMSFMDRMPLFVWS-------VLITAILLLLS-IL-WFPVLAG---AITMLLTDRNL-NTSFFDPAGGGDPTILYQHLFWFFGHPEVYILILPGFGMISHIISYYSGKKEPFGYLGMIYAMMAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKVFSWLATLHGSNIKYSPPMLWALGFIFLFTVGGLTGVVLANSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGGFIHWFPLFTGLTLNPTWLKIQFTIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIISSIGSFISLTAVMLFIFIIWEAFASKRKVLFVENTSSSIEWLQGCPPAEHSYEEPPYLKSQSNTYFIEDNWYYGVMGLCSLFNM----------------------"
        consensus = dumb_consensus([ref_msa[name]for name in ref_msa.keys()])
        t_consensus = list(consensus)

        print "Augmenting reference templates"
        # make a list of augmented msa ref seqs
        aug_msa_template_refs = dict((name, augment(msa_template_refs[name], t_consensus)) for name in name_keys)

        # Pairwise align query AA to ref AA
        print "Doing pairwise alignments..."
        i = 0
        to_go = len(name_keys)
        pairwise_refs = {}
        pairwise_queries = {}
        for name in name_keys:
            # NOTE: if you use augmented MSA templates to generage the pairwise alignement, you must also use augmented MSA
            #       templates for regapping and for delte detection.
            pairwise_refs[name], pairwise_queries[name] = globalProtAlign(
                                                            "".join(aug_msa_template_refs[name]).replace("-", ""),
                                                            translations[name][0], -2, -1)[0:2]
            i += 1
            if not i%10:
                print "%d / %d" % (i, to_go)

    # SAVE DATA
    if save_data_to_file:
        with open(savefilepath,'w') as outputfile:
            data = {"pairwise_queries":pairwise_queries,
                    "pairwise_refs":pairwise_refs,
                    "consensus": consensus
                    }
            outputfile.write(str(data))
        outputfile.close()

    """
    print "CS: " + consensus
    for name in name_keys:
        print "TR: " + str(msa_template_refs[name].seq)
        print "TA: " + "".join(aug_msa_template_refs[name])
        print "PQ: " + pairwise_queries[name]
        print "PR: " + pairwise_refs[name]
    """
    # Regap the pairwise alignments according to the original MSA
    # NOTE: if you used augmented MSA templates to generage the pairwise alignement, you must also use augmented MSA
    #       templates for regapping and delta detection.
    regappings = dict((name, nast_regap(aug_msa_template_refs[name], pairwise_refs[name], pairwise_queries[name]))
                       for name in name_keys)
    nast_queries = dict((name, regappings[name][0]) for name in name_keys)
    nast_refs = dict((name, regappings[name][1]) for name in name_keys)

    """
    print "C: " + consensus
    for name in name_keys:
        print "A: " + "".join(aug_msa_template_refs[name])
        print "R: " + nast_refs[name]
        print "Q: " + nast_queries[name]
    """
    delta_counts = {}
    # compute the deltas between the pairwise template and the MSA template
    cumulative_insertions = defaultdict(set)
    priority = 0
    for name in name_keys:
        # NOTE: if you used augmented MSA templates to generage the pairwise alignement, you must also use augmented MSA
        #       templates for regapping and delta detection.
        local_insertions, insertion_count = locate_deltas(aug_msa_template_refs[name],
                                                          nast_refs[name],
                                                          nast_queries[name],
                                                          priority)
        # Resolve the local with the set of global deltas
        update_global_deltas(local_insertions, cumulative_insertions)
        delta_counts[name] =  insertion_count
        # replace each query's deltas with lowercase so we know to replace them later (instead of inserting gaps)
        nast_queries[name] = mask_deltas(local_insertions, nast_queries[name])
        priority +=1
    print cumulative_insertions


    print "C: " + apply_deltas(cumulative_insertions, consensus, '@', True)
    # A ruler
    ruler = (' ' * 4 + "*" + ' ' * 4 + "!") * 200
    print "R: " + apply_deltas(cumulative_insertions, ruler, '@', True)

    if len(cumulative_insertions) > 0:
        for name in name_keys:
            print "T: " + apply_deltas(cumulative_insertions, msa_template_refs[name], '@')
            print "Q: " + apply_deltas(cumulative_insertions, nast_queries[name], '@')
            "\n"

    else:
        for name in name_keys:
            print "A: " + "".join(aug_msa_template_refs[name])
            print "T: " + msa_template_refs[name].seq
            print "Q: " + nast_queries[name]

    # COMPUTE REGIONS TO REALIGN
    # get_realign_segments(cumulative_insertions)
    realignment_regions = get_realignment_regions(cumulative_insertions)
    formatted_queries = [apply_deltas(cumulative_insertions, nast_queries[name], '-') for name in name_keys]

    # DO REALIGNMENT WITH MUSCLE
    i=0
    maps = []
    for (start,end) in realignment_regions:
        if end != start:
            print start, end
            maps.append(realign(formatted_queries, start, end+1, "%d_realigned.fa" % i))
            i += 1

    # Copy realigned sections back into their respective sequences
    i=0
    for (start,end) in realignment_regions:
        if end != start:
            for j in range(len(formatted_queries)):
                tmp = list(formatted_queries[j])
                tmp[start:end+1] = maps[i]["".join(tmp[start:end+1]).replace('-','').lower()].upper()
                formatted_queries[j] = tmp
            i+=1
    print "C: " + apply_deltas(cumulative_insertions, consensus, '@', True)
    for seq in formatted_queries:
        print "Q: " + "".join(seq)

    with open("%s/rslt.log" % outdir, 'w') as log:
        log.write(apply_deltas(cumulative_insertions, ruler, '@'))

    with open("%s/rslt.msa.faa" % outdir, 'w') as outputFile:
         for seq in SeqIO.parse(open(ref_faa_msa,'r'),'fasta'):
            outputFile.write(">%s\n%s\n" % (seq.id,  apply_deltas(cumulative_insertions, seq.seq, '-')))
         for query in name_keys:
            outputFile.write(">%s\n%s\n" % (query, apply_deltas(cumulative_insertions, nast_queries[query], '-')))

    with open("%s/rslt.count" % outdir, 'w') as countFile:
        for name in name_keys:
            countFile.write("%s\t%d\n" % (name, delta_counts[name]))

data_dir = "/home/greg/ARMS/src/ARMS/dev/data"
add_to_msa("%s/bold1k.fna" % data_dir,
         "%s/bold10k.fna" % data_dir,
         "%s/bold10k.faa.msa" % data_dir,
         "%s/bold10k_name_pairs.txt" % data_dir,
         "/home/greg/ARMS/src/ARMS/dev/data")



