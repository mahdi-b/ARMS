from itertools import product

from classes.ChewbaccaProgram import *
from classes.ProgramRunner import *

from annotate_OTU_table import annotateOTUtable
from classes.Helpers import *


class Annotate_OTU_Table_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"


    def execute_program(self):
        args = self.args
        self.annotate_otu_chewbacca(args.input_f, args.outdir, args.annotation, args.threads)


    def annotate_otu_chewbacca(self, input_f, outdir, annotation, threads):
        """Annotates an OTU table.

        :param input_f: Filepath to a file or folder of files to annotate.
        :param annotation: Filepath to a file or a folder of files to use as annotations.
        :param outdir: Filepath to the output directory where annotated files will be written.
        :param threads: The number of processes to use to process the input fileset.
        """
        makeDirOrdie(outdir)
        matricies = getInputFiles(input_f)
        debugPrintInputInfo(matricies, "annotated.")
        annotations = getInputFiles(annotation)
        debugPrintInputInfo(annotations, "parse.")
        inputs = product(matricies, annotations)

        pool = init_pool(min(len(matricies) * len(annotations), threads))
        printVerbose("Annotating matrix...")
        parallel(runPythonInstance, [(annotateOTUtable, matrix, annotation, "%s/%s.txt" % (outdir, "matrix"))
                                     for matrix, annotation in inputs], pool)
        printVerbose("Done Annotating.")

        cleanup_pool(pool)