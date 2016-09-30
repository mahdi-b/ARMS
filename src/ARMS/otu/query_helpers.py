from classes.ProgramRunner import *

from classes.Helpers import *


def query_vsearch(inputs, outdir, threads, aln_user_string, pool):
    """Runs a VSEARCH alignment on pairs of query/reference sequences.

    :param inputs: A list of pairs of (filepaths to) query_fastas and the refrence fastas to compare them to.
    :param outdir: Filepath to the directory where the alignment result should be written.
    :param threads: The number of threads to use per querey process.
    :param aln_user_string: An optional string of commandline parameters passed to the VSEARCH program.
    :param pool: A fully initalized multiprocessing.Pool object.
    """
    printVerbose("Aligning against reference sequences...")
    #     # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
    # --userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.ALIGN_VSEARCH,
                                                      [threads, query_fasta, ref_fasta,
                                                       "%s/%s.out" % (outdir, strip_ixes(query_fasta)),
                                                       "%s/%s.alnout" % (outdir, strip_ixes(query_fasta)),
                                                       aln_user_string],
                                                      {"exists": [query_fasta, ref_fasta], "positive": [threads]})
                                        for query_fasta, ref_fasta in inputs], pool)
    printVerbose("Done aligning.")



