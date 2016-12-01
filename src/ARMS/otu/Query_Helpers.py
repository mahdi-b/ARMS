from classes.Helpers import run_parallel, printVerbose, strip_ixes
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands

def query_vsearch(inputs, outdir, simmilarity, processes, aln_user_string, extraargstring, pool):
    """Runs a VSEARCH alignment on pairs of query/reference sequences.

    :param inputs: A list of pairs of (filepaths to) query_fastas and the refrence fastas to compare them to.
    :param outdir: Filepath to the directory where the alignment result should be written.
    :param aln_user_string: An optional string of commandline parameters passed to the VSEARCH program.
    :param simmilarity: The minimum simmilarity percentage (between reference and query sequences), \
                            as a decimal between 0 and 1), required for a positive  match.
    :param processes: The number of processes to use in the identification process.
    :param extraargstring: Advanced program parameter string.
    :param pool: A fully initalized multiprocessing.Pool object.
    """
    printVerbose("Aligning against reference sequences...")
    #     # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
    # --userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt
    run_parallel([ProgramRunner(ProgramRunnerCommands.ALIGN_VSEARCH,
                            [processes, query_fasta, ref_fasta, simmilarity, "%s/%s.out" % (outdir, strip_ixes(query_fasta)),
                             "%s/%s.alnout" % (outdir, strip_ixes(query_fasta)), aln_user_string],
                            {"exists": [query_fasta, ref_fasta], "positive": [processes]},
                            extraargstring)
              for query_fasta, ref_fasta in inputs], pool)
    printVerbose("Done aligning.")
    return