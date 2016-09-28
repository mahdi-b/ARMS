from classes.Helpers import *
from classes.ProgramRunner import *

def demux_fastx(file_id_pairs, barcodes, outdir, threads, pool):
    """Demuxes using FAST X BARCODE SPLITTER.

    :param file_id_pairs: A list of tuples of input filepaths and the corresponding shard #s. (filepath, shard#).
    :param barcodes: File path to input barcodes file.
    :param outdir: Filepath to output directory.
    :param threads: Number of processes to use to demux input fileset.
    :param pool: A fully initalized multiprocessing.Pool object.
    """
    printVerbose("Demuxing sequences...")
    parallel(runProgramRunnerInstance, [ProgramRunner(ProgramRunnerCommands.DEMUX_FASTX,
                                                      [input_, barcodes, "%s/" % outdir,
                                                       "_%d_splitOut.fastq" % id_], {"exists": [input_, barcodes]})
                                        for input_, id_ in file_id_pairs], pool)
    printVerbose("Demuxed sequences.")

    # gather output files and move them to their final destination
    output_files = enumerateDir(".", "*_splitOut_")
    bulk_move_to_dir(output_files, outdir)

    # Grab all the auxillary files
    aux_files = getInputFiles(outdir, "unmatched_*")
    # make aux dir for extraneous files and move them there
    bulk_move_to_dir(aux_files, makeAuxDir(outdir))