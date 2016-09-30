from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *

from classes.Helpers import *


class Preclean_Program_Bayeshammer(ChewbaccaProgram):
    """Uses bayeshammer (Spades) to fix sequencing errors via kmer clustering and probabilistic substitution.
    """
    name = "bayeshammer"


    def execute_program(self):
        args = self.args
        self.preclean_bayeshammer(args.input_f, args.input_r, args.outdir, args.threads, args.bayesthreads)


    def preclean_bayeshammer(self, input_f, input_r, outdir, threads, bayesthreads):
        """Assembles reads from two (left and right) fastq files/directories.

        :param input_f: File path to file or folder of left reads to clean.
        :param input_r: File path to file or folder of right reads to clean.
        :param outdir: Filepath to output directory.
        :param program: The program to use.  Current options are "bayeshammer".
        :param threads: Number of processes to use to process the input files.
        :param bayesthreads: The number of threads per process to use.
        """
        makeDirOrdie(outdir)
        # Collect input files, and validate that they match
        inputs = validate_paired_fastq_reads(input_f, input_r)
        pool = init_pool(min(len(inputs), threads))
        printVerbose("\tPrecleaning reads with Spades-Baye's Hammer...")
        debugPrintInputInfo(inputs, "preclean/fix.")

        parallel(runProgramRunnerInstance,
                 [ProgramRunner(ProgramRunnerCommands.PRECLEAN_SPADES,
                                [forwards, reverse, outdir, bayesthreads],
                                {"exists": [forwards, reverse], "positive": [bayesthreads]})
                    for forwards, reverse in inputs], pool)
        printVerbose("Done cleaning reads.")

        # Grab all the auxillary files (everything not containing ".assembled."
        #aux_files = getInputFiles(outdir, "*", "*.assembled.*", ignore_empty_files=False)
        # make aux dir for extraneous files and move them there
        #bulk_move_to_dir(aux_files, makeAuxDir(outdir))

        # Select output files
        cleaned_reads = getInputFiles("%s/corrected" % (outdir), "*.gz")
        bulk_move_to_dir(cleaned_reads, outdir)


        # Gather aux files
        aux_dir = makeAuxDir(outdir)
        aux_files = getInputFiles(outdir, "*", "*.gz")
        aux_files += getInputFiles(outdir, "*unpaired*")
        bulk_move_to_dir(aux_files, aux_dir)
        cleanup_pool(pool)