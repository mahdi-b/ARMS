from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *

from classes.Helpers import *
from partition import splitK


class Partition_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"


    def execute_program(self):
        args = self.args
        self.partition_chewbacca(args.input_f, args.outdir, args.threads, args.chunksize, args.filetype)


    def partition_chewbacca(self, input_f, outdir, threads, chunksize, filetype):
        """Partition a fasta/fastq file into chunks of user-defined size.

        :param input_f: Filepath to a file or folder of files to partition.
        :param outdir: The directory to write split files to.
        :param threads: The number of processes to use to partition the input fileset.
        :param chunksize: The number of sequences per file.
        :param filetype: Either 'fasta' or 'fastq'.
        """
        # def splitK(inputFasta, prefix, nbSeqsPerFile, filetype):
        makeDirOrdie(outdir)
        # Gather input files
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "partitioned")
        pool = init_pool(min(len(inputs), threads))
        printVerbose("Partitioning Files...")
        parallel(runPythonInstance,
                 [(splitK, input_, "%s/%s" % (outdir, strip_ixes(input_)), chunksize, filetype)
                  for input_ in inputs], pool)
        printVerbose("Done partitioning files.")
        cleanup_pool(pool)