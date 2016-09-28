from classes.Helpers import *
from classes.ProgramRunner import *


def preclean_bayeshammer(inputs, outdir, threads, pool, spadesthreads):
    """Assembles reads from two (left and right) fastq files/directories.  For a set of k forward read files, and k
        reverse read files, return k assembled files.  Matching forward and reverse files should be identically named,
        except for a <forward>/<reverse> suffix that indicates the read orientation.  Two suffix pairs are supported:
        '_forwards' and '_reverse',
        and
        '_R1' and 'R2'
        Choose ONE suffix style and stick to it.
        e.g. Sample_100_forwards.fq and Sample_100_reverse.fq will be assembled into Sample_100_assembled.fq.
          Alternatively, Sample_100_R1.fq and Sample_100_R2.fq will be assembled into Sample_100_assembled.fq.
          You can provide as many pairs of files as you wish as long as they follow exactly on of the above naming
          conventions.  If a 'name' parameter is provided, it will be used as a suffix for all assembled sequence files.
    :param args: An argparse object with the following parameters:
                    name            Textual ID for the data set.
                    input_f         Forward Fastq Reads file or directory.
                    input_r         Reverse Fastq Reads file or directory.
                    threads         The number of threads to use durring assembly.
                    outdir          Directory where outputs will be saved.
    """

    printVerbose("\tPrecleaning reads with Spades-Baye's Hammer...")
    debugPrintInputInfo(inputs, "preclean/fix.")

    parallel(runProgramRunnerInstance,
             [ProgramRunner(ProgramRunnerCommands.PRECLEAN_SPADES,
                            [forwards, reverse, outdir, spadesthreads],
                            {"exists": [forwards, reverse], "positive": [spadesthreads]})
                for forwards, reverse in inputs], pool)
    printVerbose("Done cleaning reads.")

    # Grab all the auxillary files (everything not containing ".assembled."
    #aux_files = getInputFiles(outdir, "*", "*.assembled.*", ignore_empty_files=False)
    # make aux dir for extraneous files and move them there
    #bulk_move_to_dir(aux_files, makeAuxDir(outdir))

    # Select output files
    cleaned_reads = getInputFiles("%s/corrected" % (outdir), "*.gz")
    bulk_move_to_dir(cleaned_reads, outdir)

    aux_dir = makeAuxDir(outdir)
    aux_files = getInputFiles(outdir, "*", "*.gz")
    aux_files += getInputFiles(outdir, "*unpaired*")
    bulk_move_to_dir(aux_files, aux_dir)
    # Gather aux files