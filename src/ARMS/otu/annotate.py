from classes.Helpers import *
from classes.ProgramRunner import *
from annotateOTUtable import annotateOTUtable
from itertools import product

def annotate_main(input_f, annotation, outdir, threads, program, aux_params):
    """Annotates an OTU table.

    :param input_f: Filepath to a file or folder of files to annotate.
    :param annotation: Filepath to a file or a folder of files to use as annotations.
    :param outdir: Filepath to the output directory where annotated files will be written.
    :param threads: The number of processes to use to process the input fileset.
    :param program: The Program to use to annotate the OTU tables.  Choices are ["chewbacca"].  Default: "chewbacca"
    :param aux_params: A dictionary of program-specific named-parameters.
    """
    makeDirOrdie(outdir)
    matricies = getInputFiles(input_f)
    debugPrintInputInfo(matricies, "annotated.")
    annotations = getInputFiles(annotation)
    debugPrintInputInfo(annotations, "parse.")
    inputs = product(matricies, annotations)

    pool = init_pool(min(len(matricies) * len(annotations), threads))
    if program == "chewbacca":
        annotate_chewbacca(inputs, outdir, pool)
    cleanup_pool(pool)
    
    
def annotate_chewbacca(inputs, outdir, pool):
    """Given a list of paired OTU tables and annotation files, use the annotation files to annotate the associated
        OTU table.

    :param inputs: A list of pairs of OTU tables and annotation files for those tables.
    :param outdir: The output directory to write annotated OTU tables to.
    :param pool:  A fully-initalized multiprocessing.Pool object.
    """
    printVerbose("Annotating matrix...")
    parallel(runPythonInstance, [(annotateOTUtable, matrix, annotation, "%s/%s.txt" % (outdir, "matrix"))
                                 for matrix, annotation in inputs], pool)
    printVerbose("Done Annotating.")
