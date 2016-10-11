from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import getInputFiles, debugPrintInputInfo, init_pool, run_parallel, printVerbose, strip_ixes, \
    cleanup_pool, bulk_move_to_dir, makeAuxDir
from classes.PythonRunner import PythonRunner
from itertools import product
from queryVSearchoutForTaxa import parseVSearchOutputAgainstFasta
from Query_Helpers import query_vsearch


class Query_OTU_Fasta_Program_Vsearch(ChewbaccaProgram):
    name = "vsearch"

    def execute_program(self):
        args = self.args
        self.query_fasta_vsearch(args.input_f, args.referencefasta, args.taxinfo, args.outdir, args.processes,
                                 args.simmilarity, args.coverage, args.extraargstring)

    def query_fasta_vsearch(self, input_f, referencefasta, taxinfo, outdir, processes, simmilarity, coverage,
                            extraargstring):
        """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

        :param input_f:  Filepath to a file or folder of files to identify.
        :param outdir: Filepath to the output directory.
        :param referencefasta: Filepath to a file or folder of files to use as a reference.
        :param taxinfo:  Filepath to a file containing taxonomic info correlated with the referencefasta.
        :param simmilarity: The % simmilarity between a query and reference sequence required for positive
                                identification.
        :param coverage: The % coverage of matching regions between a query and reference sequence required for positive
                            identification.
        :param processes: The number of processes to use in the identification process.
        :param extraargstring: Advanced program parameter string.
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
            print "Error: The number of reference fastas and taxonomic mapping files is not the same.  There must be \
                    one taxonomic mapping file for each reference fasta."
            return
        ref_data_pairs = zip(ref_fastas, tax_info_files)
        inputs = [x for x in product(query_fastas, ref_fastas)]
        aln_user_string = ""
        pool = init_pool(min(len(inputs), processes))

        # VSEARCH ALIGNMENT
        query_vsearch(inputs, outdir, processes, aln_user_string, extraargstring, pool)

        printVerbose("Parsing output...")
        # Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in
        # parsed_BIOCODE.out.  Parameters can be changed and this command can be rerun as many times as necessary
        #
        # parseVSearchOutputAgainstFasta(vsearch_outfile, taxInfo, output_file, min_simmilarity, min_coverage):
        inputs = [x for x in product(query_fastas, ref_data_pairs)]
        debugPrintInputInfo(inputs, "queryied against paired refereces.")
        run_parallel([PythonRunner(parseVSearchOutputAgainstFasta,
                                   ["%s/%s.out" % (outdir, strip_ixes(query)), tax_info,
                                    "%s/%s.tax" % (outdir, strip_ixes(query)), simmilarity, coverage],
                                   {"exists": [query, ref_fasta, tax_info]})
                      for query, (ref_fasta, tax_info) in inputs], pool)
        printVerbose("\nDone parsing...")

        # Gather and move auxillary files
        aux_files = getInputFiles(outdir, "*", "*.tax", ignore_empty_files=False)
        bulk_move_to_dir(aux_files, makeAuxDir(outdir))

        cleanup_pool(pool)
