from itertools import product

from classes.ProgramRunner import *
from enum import Enum
from parse.parseVSearchoutForTaxa import parseVSearchOutputAgainstFasta, parseVSearchOutputAgainstNCBI

from classes.Helpers import *


class QUERY_TYPE(Enum):
    FASTA = "FASTA"
    DATABASE = "DATABASE"

def query_main(input_f, outdir, threads, program, ref, aux_params):
    """Performs closed-reference OTU picking (alignment against known reference sequences).  Switchboard function.

    :param input_f:  Filepath to a file or folder of files to identify.
    :param outdir: Filepath to the output directory.
    :param threads: The number of threads to use in the identification process.
    :param program: The identification program to use.  Choices are ["vsearch"].  Default: "vsearch".
    :param ref: The reference type.  Either "fasta" or "NCBIDB".
    :param aux_params: A dictionary of program-specific named-parameters.
    """
    makeDirOrdie(outdir)
    if program == "vsearch":
        if ref == QUERY_TYPE.DATABASE:
            query_ncbi_vsearch(input_f, outdir, threads)
        elif ref == QUERY_TYPE.FASTA:
            query_fasta_vsearch(input_f, aux_params["referencefasta"], aux_params["taxinfo"], outdir, threads,
                            aux_params["simmilarity"], aux_params["coverage"])


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


def query_fasta_vsearch(input_f, referencefasta, taxinfo, outdir, threads, simmilarity, coverage):
    """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

    :param input_f:  Filepath to a file or folder of files to identify.
    :param outdir: Filepath to the output directory.
    :param threads: The number of threads to use in the identification process.
    :param referencefasta: Filepath to a file or folder of files to use as a reference.
    :param taxinfo:  Filepath to a file containing taxonomic info correlated with the referencefasta.
    :param simmilarity: The % simmilarity between a query and reference sequence required for positive identification.
    :param coverage: The % coverage of matching regions between a query and reference sequence required for positive
                        identification.
    """
    # vsearch --usearch_global %s seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
    #       --userfields query+target+id+alnlen+qcov --userout %sout --alnout %s alnout.txt

    # expecting a fasta to annotate
    query_fastas = getInputFiles(input_f)
    debugPrintInputInfo(query_fastas, "queried for identification.")
    ref_fastas = getInputFiles(referencefasta)
    debugPrintInputInfo(ref_fastas, "referenced for sequence identification.")
    tax_info_files = getInputFiles(taxinfo)
    debugPrintInputInfo(tax_info_files, "referenced for taxanomic names.")

    # make sure the number of reference fasta files is the same as the number of tax_info files
    if len(tax_info_files) != len(ref_fastas):
        print "Error: The number of reference fastas and taxonomic mapping files is not the same.  There must be one \
                taxonomic mapping file for each reference fasta."
        return
    ref_data_pairs = zip(ref_fastas, tax_info_files)
    inputs = [x for x in product(query_fastas, ref_fastas)]
    aln_user_string = ""
    pool = init_pool(min(len(inputs), threads))

    # VSEARCH ALIGNMENT
    query_vsearch(inputs, outdir, threads, aln_user_string, pool)

    printVerbose("Parsing output...")
    # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
    # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
    #
    # parseVSearchOutputAgainstFasta(vsearch_outfile, taxInfo, output_file, min_simmilarity, min_coverage):
    inputs = [x for x in product(query_fastas, ref_data_pairs)]
    parallel(runPythonInstance, [(parseVSearchOutputAgainstFasta, "%s/%s.out" % (outdir, strip_ixes(query)),
                                  tax_info, "%s/%s_result.out" % (outdir, strip_ixes(query)),
                                  simmilarity, coverage)
                                 for query, (ref_fasta, tax_info) in inputs], pool)
    printVerbose("\nDone parsing...")
    # Gather and move auxillary files
    aux_files = getInputFiles(outdir, "*", "*_result.out")
    bulk_move_to_dir(aux_files, makeAuxDir(outdir))

    cleanup_pool(pool)


def query_ncbi_vsearch(input_f, outdir, threads):
    """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

    :param input_f:  Filepath to a file or folder of files to identify.
    :param outdir: Filepath to the output directory.
    :param threads: The number of threads to use in the identification process.
    """
    coi_fasta = os.path.expanduser("~/ARMS/refs/COI.fasta")
    ncbi_db_string = os.path.expanduser("~/ARMS/refs/ncbi.db")
    aln_user_string = "--userfields query+target+id+alnlen+qcov"

    query_fastas = getInputFiles(input_f)
    debugPrintInputInfo(query_fastas, "queried against NCBI.")
    inputs = [x for x in product(query_fastas, [coi_fasta])]
    pool = init_pool(min(len(query_fastas), threads))

    # VSEARCH ALIGNMENT
    query_vsearch(inputs, outdir, threads, aln_user_string, pool)

    printVerbose("Parsing output...")
    # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
    # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
    #
    # parseVSearchOutputAgainstNCBI(vsearch_out, ncbi_db, min_coverage, min_similarity)> parsed_nt.out
    parallel(runPythonInstance,
             [(parseVSearchOutputAgainstNCBI, "%s/%s.out" % (outdir, strip_ixes(query)), ncbi_db_string,
               "%s/%s_result.out" % (outdir, strip_ixes(query)),
               97, 85) for query in query_fastas], pool)
    printVerbose("Done processing.")

    # Gather and move auxillary files
    aux_dir = makeAuxDir(outdir)
    aux_files = getInputFiles(outdir, "*", "*_result.out")
    bulk_move_to_dir(aux_files, aux_dir)

    cleanup_pool(pool)
