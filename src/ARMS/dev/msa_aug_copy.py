from Bio import SeqIO
from collections import defaultdict
from nast import nast_regap, locate_deltas, apply_deltas, mask_deltas, make_faa_gc_lookup, \
    get_best_hits_from_vsearch, translate_to_prot, update_global_deltas, get_realign_segments, get_region, augment, realign
from utils import globalProtAlign

import os

#query.fna", bold10k.fna", bold10k.faa.msa, bold10k_name_pairs.txt
def add_to_msa(input_fna, ref_fna, ref_faa_msa, name_map, outdir):

    # read MSA file to dict
    print "Reading MSA file"
    ref_msa = SeqIO.to_dict(SeqIO.parse(open(ref_faa_msa, 'r'), 'fasta'))

    # store the gc for each ref seq
    print "Looking up GCs and names"
    gc_lookup_map, fna_faa_map = make_faa_gc_lookup(name_map)

    # Find closest neighbors for each input sequence
    print "Finding best pairs"
    best_hits = get_best_hits_from_vsearch(input_fna, ref_fna, outdir)
    if len(best_hits) == 0:
        print "\n\n\nERROR: no matching IDs found in reference DB."
        exit()
    # Translate NA to AA
    print "Translating"
    translations = translate_to_prot(input_fna, best_hits, gc_lookup_map)





    # Regap pairwise alignment to include MSA template gaps
    msa_templates = []
    nast_refs = []
    nast_queries = []
    out_queries = {}
    pairwise_queries = []
    pairwise_refs = []
    cumulative_insertions = defaultdict(set)
    priority = 0
    i = 0
    consensus = "-------------------------------------G---MMVMFADRWLFSTNHKDIGTLYFIFGAWIAGM--IVGTSLSL------------------LIRAELGQPPG-SHLLNI----GDNDQ--S--FL--------------IYNVIVTAHAFIMIFF--MVMPIMIGGFG--L-------NWLVPLMLNTSFFDPAGGAPDMAFPRMNNMS---F---K--FW--LLGHPPSLTLLLSSS-----M-----I-VENGRKK-F--KC-----TA-------KKEAFGTGWTVYP---PLSSNILGFVVWAHSGASVD----------------DLAI------FS----------------------LHLAGISSILG---------------------------AINFI--------------TTIINMSR--IISNGMSFMDRMPLFVWS-------VLITAILLLLS-IL-WFPVLAG---AITMLLTDRNL-NTSFFDPAGGGDPTILYQHLFWFFGHPEVYILILPGFGMISHIISYYSGKKEPFGYLGMIYAMMAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKVFSWLATLHGSNIKYSPPMLWALGFIFLFTVGGLTGVVLANSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGGFIHWFPLFTGLTLNPTWLKIQFTIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIISSIGSFISLTAVMLFIFIIWEAFASKRKVLFVENTSSSIEWLQGCPPAEHSYEEPPYLKSQSNTYFIEDNWYYGVMGLCSLFNM----------------------"
    savefile = "%s.save" % input_fna
    loaded=False
    if os.path.isfile(savefile) and False:
        lines = "".join([line.rstrip() for line in open(savefile,'r')])
        loaded = True
        data = eval(lines)
    print consensus

    for name in translations.keys():
        msa_template_ref = ref_msa[fna_faa_map[best_hits[name][1]]]
        if loaded:
            pairwise_ref = data["pairwise_refs"][i]
            pairwise_query = data["pairiwse_queries"][i]
        else:
            # Pairwise align query AA to ref AA
            t_consensus = list(consensus)
            t_msa_ref = msa_template_ref
            t_msa_ref = augment(t_msa_ref, t_consensus)
            pairwise_ref, pairwise_query = globalProtAlign("".join(t_msa_ref).replace("-",""), translations[name][0])[0:2]

        # Replace template gaps with consensus data
        augment(pairwise_ref, pairwise_query)

        pairwise_queries.append(str(pairwise_query))
        pairwise_refs.append(str(pairwise_ref))
        print "S: " + str(msa_template_ref.seq).replace("-","")
        print "A: " + pairwise_ref
        print "Q: " + pairwise_query

        nast_query, nast_ref = nast_regap(msa_template_ref.seq, pairwise_ref, pairwise_query)

        i += 1


        # trim msa ref and nast ref strings to same size
        nast_ref_len = len(nast_ref)
        msa_temp_len = len(msa_template_ref)
        """
        if nast_ref_len > msa_temp_len:
            nast_query = nast_query[nast_ref_len - msa_temp_len:]
            nast_ref = nast_ref[nast_ref_len - msa_temp_len:]
        if msa_temp_len > nast_ref_len:
            nast_query = '-' * (msa_temp_len - nast_ref_len) + nast_query
            nast_ref = '-' * (msa_temp_len - nast_ref_len) + nast_ref
        """
        # compute the deltas between the pairwise template and the MSA template
        local_insertions, total_insertions = locate_deltas(msa_template_ref, nast_ref, nast_query, priority)

        print "M: " + msa_template_ref
        print "R: " + nast_ref
        print "Q: " + nast_query

        if total_insertions < 10:
            update_global_deltas(local_insertions, cumulative_insertions)
            # replace each query's deltas with lowercase so we know to replace them later (instead of inserting gaps)
            nast_query = mask_deltas(local_insertions, nast_query)

            priority +=1
            # nast_refs.append(nast_ref)
            nast_queries.append(nast_query)
            msa_templates.append(msa_template_ref)
            out_queries[name] = nast_query
        i+=1
        if i%10 ==0: print "%d"%i

    if False:
        with open(savefile,'w') as outputfile:
            data = {"pairiwse_queries":pairwise_queries,
                    "pairwise_refs":pairwise_refs
                    }
            outputfile.write(str(data))
        outputfile.close()

    print cumulative_insertions

    # A ruler
    ruler = (' ' * 4 + "*" + ' ' * 4 + "!") * 200
    get_realign_segments(cumulative_insertions)

    print "R: " + apply_deltas(cumulative_insertions, ruler, '@', True)

    if len(cumulative_insertions) > 0:

        for i in range(len(msa_templates)):
            print "T: " + apply_deltas(cumulative_insertions, msa_templates[i], '@')
            print "Q: " + apply_deltas(cumulative_insertions, nast_queries[i], '@')
            "\n"
        # insert back into msa, and update msa

    else:
        for i in range(len(msa_templates)):
            print msa_templates[i]
            print nast_queries[i]
            print "\n"


    segments = get_region(cumulative_insertions, ruler)
    formatted_queries = [apply_deltas(cumulative_insertions, seq, '-') for seq in nast_queries]

    # DO REALIGNMENT WITH MUSCLE
    i=0
    maps = []
    segments = [(0,150)]
    for (start,end) in segments:
        if end != start:
            maps.append(realign(formatted_queries, start, end+1, "%d_realigned.fa" % i))
            i += 1
    i=0
    for (start,end) in segments:
        if end != start:
            for j in range(len(formatted_queries)):
                tmp = list(formatted_queries[j])
                tmp[start:end+1] = maps[i]["".join(tmp[start:end+1]).lower()].upper()
                formatted_queries[j] = tmp
            i+=1
    for seq in formatted_queries:
        print "".join(seq)

    with open("%s/rslt.log" % outdir, 'w') as log:
        log.write(apply_deltas(cumulative_insertions, ruler, '@'))

    with open("%s/rslt.faa.msa" % outdir, 'w') as outputFile:
        for query in out_queries.keys():
            outputFile.write(">%s\n%s\n" % (query, apply_deltas(cumulative_insertions, out_queries[query], '@')))

# get a list of @@ segments
# write to file: each @@ segment.
# unique each @@ segment
# call muscle on each segment
# do map of old segments to new segments
# reconstruct each sequence
    # lookup each @@ segment in the dict to get the new segment
# write to out file
data_dir = "/home/greg/ARMS/src/ARMS/dev/data"
add_to_msa("%s/bold50.fna" % data_dir,
         "%s/bold10k.fna" % data_dir,
         "%s/bold10k.faa.msa" % data_dir,
         "%s/bold10k_name_pairs.txt" % data_dir,
         "/home/greg/ARMS/src/ARMS/dev/data")



