from itertools import product

from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *
from parse.parseVSearchoutForTaxa import parseVSearchOutputAgainstNCBI

from classes.Helpers import *
from query_helpers import query_vsearch


class Query_OTU_DB_Program_Vsearch(ChewbaccaProgram):
    name = "vsearch"


    def execute_program(self):
        args = self.args
        self.query_fasta_db_vsearch(args.input_f, args.outdir, args.referencefasta, args.db,  args.threads)



    def query_fasta_db_vsearch(self, input_f, outdir, ref_fasta, ref_db, threads):
        """Compare reference sequences to the fasta-formatted query sequences, using global pairwise alignment.

        :param input_f:  Filepath to a file or folder of files to identify.
        :param outdir: Filepath to the output directory.
        :param threads: The number of threads to use in the identification process.
        """
        makeDirOrdie(outdir)
        aln_user_string = "--userfields query+target+id+alnlen+qcov"
        #coi_fasta = os.path.expanduser("~/ARMS/refs/COI.fasta")
        #ncbi_db_string = os.path.expanduser("~/ARMS/refs/ncbi.db")
        coi_fasta = ref_fasta
        ncbi_db_string = ref_db

        query_fastas = getInputFiles(input_f)
        debugPrintInputInfo(query_fastas, "queried against the DB.")
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
