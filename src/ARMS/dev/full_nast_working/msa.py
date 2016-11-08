from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from classes.Helpers import init_pool, strip_ixes, printVerbose
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from itertools import product
from dev.utils import globalProtAlign
import subprocess
"""
How to get an MSA for all new seqs against the old seqs?

    Input:
    ------
        ref unaligned fnas (and some annotation for gc)
        ref AA MSA alignment
        query seq

    Process:
    --------
    Step1: ID best match:
        vsearch to align query against unaligned fnas
        get best vsearch reslut,
    Step2: Translate NA to AA
        translate query to AA using annotated gc from that match
    Step3: Align query AA to ref AA
        align query AA to ref AA, using pairwise alignment; only insert gaps into query seq (don't mess up the ref MSA)

    Output:
    -------
        MSA of queries

Result: All querys should then align to the ref MSA by transitivity!
"""




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
    # Search for good hits
    inputs = [(input_fna, ref_fna)]
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
    # for each guy in best_hits,
        # lookup its best hit
        # lookup the gc of its best hit
        # translate into all orf using gc
        # store the one without stops
    # return translated prots

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


def make_msa(input_fna, ref_fna, ref_faa_msa, name_map, outdir):
    # read MSA file to dict
    print "Reading MSA file"
    ref_msa = SeqIO.to_dict(SeqIO.parse(open(ref_faa_msa, 'r'), 'fasta'))

    print "Reading query file"
    queries = SeqIO.to_dict(SeqIO.parse(open(input_fna, 'r'), 'fasta'))

    # store the gc for each ref seq
    print "Looking up GCs and names"
    gc_lookup_map, fna_faa_map = make_faa_gc_lookup(name_map)

    # Find closest neighbors for each input sequence
    print "Finding best pairs"
    best_hits = get_best_hits_from_vsearch(input_fna, ref_fna, outdir)

    # Translate NA to AA
    print "Translating"
    translations = translate_to_prot(input_fna, best_hits, gc_lookup_map)

    # Align query AA to ref AA
    print "Doing pairwise alignments"
    alignments = dict([(name, globalProtAlign(ref_msa[fna_faa_map[best_hits[name][1]]].seq.ungap('-'), translations[name][0])) for
                  name in translations.keys()])

    f_str = "Template MSA: \n=============\n%s\n\nQuery:\n======\n%s\n\nPairwise Ref:\n=============\n%s\n\nPairwise \
                Query:\n===============\n%s\n\n\n\n"
    cmd = "./dev/a.out %s %s %s %s"
    output = []
    for name in alignments.keys():
        # ref_msa, aligned query
        template_ref = ref_msa[fna_faa_map[best_hits[name][1]]].seq
        query = queries[name].seq
        pairwise_ref = alignments[name][0]
        pairwise_query = alignments[name][1]
        #print  f_str % (template_ref, query, pairwise_ref, pairwise_query)
        #"cand_ref cand_aln temp_ref temp_aln""
        cmd_str = cmd % (query, pairwise_query, template_ref, pairwise_ref)
        print "\n\n\n\nOUTPUT:\n"
        subprocess.check_call(cmd_str, shell=True)




data_dir = "/home/greg/ARMS/src/ARMS/dev/data"
make_msa("%s/query.fna"% data_dir,
         "%s/bold10k.fna" % data_dir,
         "%s/bold10k.faa.msa" % data_dir,
         "%s/bold10k_name_pairs.txt" % data_dir,
         "/home/greg/ARMS/testARMS/out1_")



