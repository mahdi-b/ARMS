from classes.Helpers import *
from classes.ProgramRunner import *

def clean_quality_trimmomatic(inputs, outdir, window_size, quality, min_len, pool):
    """Uses a sliding window to identify and trim away areas of low quality.

    :param inputs: A list of input files to clean.
    :param outdir: Filepath to the output directory.
    :param window_size: Width of the sliding window. (Number of consecutive base-pairs to average for quality analysis).
    :param quality: Minimum quality allowed.  Sections with lower average quality than this will be dropped.
    :param min_len: Minimum allowed length for TRIMMED sequences.  (i.e. if a sequence is too short after trimming,
                    its dropped.)
    :param pool: A fully-initalized multiprocessing.Pool object.
    """
    # "trimomatic":       "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
    # -%phred %input %output SLIDINGWINDOW:%windowsize:%minAvgQuality MINLEN:%minLen"

    printVerbose("Cleaning sequences with Trimmomatic...")
    parallel(runProgramRunnerInstance,
             [ProgramRunner(ProgramRunnerCommands.CLEAN_TRIMMOMATIC,
                            [input_, "%s/%s_cleaned.fastq" % (outdir, strip_ixes(input_)), window_size, quality,
                             min_len],
                            {"exists": [outdir, input_], "positive": [window_size, quality, min_len]})
              for input_ in inputs], pool)
    printVerbose("Done cleaning sequences.")