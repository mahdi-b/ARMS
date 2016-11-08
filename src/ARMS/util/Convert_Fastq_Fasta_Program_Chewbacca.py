from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.PythonRunner import PythonRunner
from convert_Fastq_To_Fasta import translateFastqToFasta

from classes.Helpers import getInputFiles, debugPrintInputInfo, printVerbose, init_pool, run_parallel, getFileName, \
                                cleanup_pool


class Convert_Fastq_Fasta_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"

    def execute_program(self):
        args = self.args
        self.convert_chewbacca(args.input_f, args.outdir, args.processes)

    def convert_chewbacca(self, input_f, outdir, proceses):
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "convert to fasta.")
        printVerbose("Converting to fasta...")
        pool = init_pool(min(len(inputs), proceses))
        run_parallel([PythonRunner(translateFastqToFasta,
                                   [input_, "%s/%s.fasta" % (outdir, getFileName(input_))],
                                   {"exists": input_})
                                    for input_ in inputs], pool)
        printVerbose("Done converting.")
        cleanup_pool(pool)