from classes.Helpers import *
from classes.ProgramRunner import *
from classes.ChewbaccaProgram import *

class Assemble_Program_Pear(ChewbaccaProgram):
    name = "pear"


    def execute_program(self):
        args = self.args
        self.assemble_pear(args.input_f, args.input_r, args.outdir, args.name,  args.threads, args.pearthreads)


    def assemble_pear(self, input_f, input_r, outdir, name, threads, pearthreads):
        """Uses PEAR to assemble paired F/R read files in parallel.

        :param input_f: File path to forward Fastq Reads file or folder.
        :param input_r: File path to reverse Fastq Reads file or folder.
        :param outdir: File path to the output directory.
        :param name: File prefix for the assembled reads.
        :param threads: The number of processes to use.
        :param pearthreads: The number of threads per process to use.
        """
        # "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s -m %d"
        makeDirOrdie(outdir)
        inputs = validate_paired_fastq_reads(input_f, input_r)
        pool = init_pool(min(len(inputs), threads))
        printVerbose("\tAssembling reads with pear")
        debugPrintInputInfo(inputs, "assemble")
        parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.ASSEMBLE_PEAR, [forwards, reverse,
                                                                                                "%s/%s_%s" % (
                                                                                                    outdir, name,
                                                                                                    getFileName(forwards)),
                                                                                                pearthreads],
                                                          {"exists": [forwards, reverse], "positive": [pearthreads]})
                                            for forwards, reverse in inputs], pool)

        printVerbose("Done assembling sequences...")
        # Grab all the auxillary files (everything not containing ".assembled."
        aux_files = getInputFiles(outdir, "*", "*.assembled.*", ignore_empty_files=False)
        # make aux dir for extraneous files and move them there
        bulk_move_to_dir(aux_files, makeAuxDir(outdir))
        cleanup_pool(pool)