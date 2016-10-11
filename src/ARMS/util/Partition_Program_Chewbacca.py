from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import getInputFiles, debugPrintInputInfo, init_pool, run_parallel, printVerbose, strip_ixes, \
    cleanup_pool
from classes.PythonRunner import PythonRunner
from partition import splitK


class Partition_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"


    def execute_program(self):
        args = self.args
        self.partition_chewbacca(args.input_f, args.outdir, args.processes, args.chunksize, args.filetype)


    def partition_chewbacca(self, input_f, outdir, processes, chunksize, filetype):
        """Partition a fasta/fastq file into chunks of user-defined size.

        :param input_f: Filepath to a file or folder of files to partition.
        :param outdir: The directory to write split files to.
        :param processes: The number of processes to use to partition the input fileset.
        :param chunksize: The number of sequences per file.
        :param filetype: Either 'fasta' or 'fastq'.
        """
        # Gather input files
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "partitioned")
        pool = init_pool(min(len(inputs), processes))
        printVerbose("Partitioning Files...")
        run_parallel([PythonRunner(splitK,
                                   [input_, "%s/%s" % (outdir, strip_ixes(input_)), chunksize, filetype],
                                   {"exists": [input_]})
                  for input_ in inputs], pool)
        printVerbose("Done partitioning files.")
        cleanup_pool(pool)