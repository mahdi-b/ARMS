from annotate_OTU_table import annotateOTUtable
from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import getInputFiles, debugPrintInputInfo, init_pool, run_parallel, printVerbose, copy_file, \
    cleanup_pool
from classes.PythonRunner import PythonRunner
from itertools import product


class Annotate_OTU_Table_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"

    def execute_program(self):
        args = self.args
        self.annotate_otu_chewbacca(args.input_f, args.outdir, args.annotation, args.processes)

    def annotate_otu_chewbacca(self, input_f, outdir, annotation, processes):
        """Annotates an OTU table.

        :param input_f: Filepath to a file or folder of files to annotate.
        :param annotation: Filepath to a file or a folder of files to use as annotations.
        :param outdir: Filepath to the output directory where annotated files will be written.
        :param processes: The maximum number of processes to use.
        """
        matricies = getInputFiles(input_f)
        debugPrintInputInfo(matricies, "annotated.")
        annotations = getInputFiles(annotation)
        # if all the annotations files are empty, just copy over files.
        if len(annotations) == 0 and len(getInputFiles(annotation, ignore_empty_files=False)) > 0:
            pool = init_pool(min(len(matricies), processes))
            print "**WARNING**: Annotation File is empty.  Skipping annotation and copying old OTU tables to output \
                    directory.\n"
            run_parallel([PythonRunner(copy_file, [matrix, outdir],
                                       {"exists": [matrix]}) for matrix in matricies], pool)
        else:
            pool = init_pool(min(len(matricies) * len(annotations), processes))
            debugPrintInputInfo(annotations, "parsed.")
            inputs = product(matricies, annotations)

            printVerbose("Annotating matrix...")
            run_parallel([PythonRunner(annotateOTUtable, [matrix, annotation, "%s/%s.txt" % (outdir, "matrix")],
                                       {"exists": [matrix, annotation]})
                          for matrix, annotation in inputs], pool)
            printVerbose("Done Annotating.")

        cleanup_pool(pool)
