from classes.Helpers import *
from classes.ProgramRunner import *
from classes.ChewbaccaProgram import *

class Clean_Quality_Program_Trimmomatic(ChewbaccaProgram):
    name = "trimmomatic"


    def execute_program(self):
        args = self.args
        self.clean_quality_trimmomatic(args.input_f, args.outdir, args.windowSize, args.quality, args.minlen,
                                       args.threads)


    def clean_quality_trimmomatic(self, input_f, outdir, window_size, quality, min_len, threads):
        """Uses a sliding window to identify and trim away areas of low quality.

        :param input_f: Filepath to input file or folder.
        :param outdir: Filepath to the output directory.
        :param window_size: Width of the sliding window. (Number of consecutive base-pairs to average for quality analysis).
        :param quality: Minimum quality allowed.  Sections with lower average quality than this will be dropped.
        :param min_len: Minimum allowed length for TRIMMED sequences.  (i.e. if a sequence is too short after trimming,
                        its dropped.)
        :param threads: Number of processes to use to clean the input fileset.
        """
        # "trimomatic":       "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
        # -%phred %input %output SLIDINGWINDOW:%windowsize:%minAvgQuality MINLEN:%minLen"
        makeDirOrdie(outdir)
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "clean")
        pool = init_pool(min(len(inputs), threads))

        printVerbose("Cleaning sequences with Trimmomatic...")
        parallel(runProgramRunnerInstance,
                 [ProgramRunner(ProgramRunnerCommands.CLEAN_TRIMMOMATIC,
                                [input_, "%s/%s_cleaned.fastq" % (outdir, strip_ixes(input_)), window_size, quality,
                                 min_len],
                                {"exists": [outdir, input_], "positive": [window_size, quality, min_len]})
                  for input_ in inputs], pool)
        printVerbose("Done cleaning sequences.")
        cleanup_pool(pool)
