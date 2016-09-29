from classes.ChewbaccaProgram import *
from classes.Helpers import *
from classes.ProgramRunner import *


class Build_OTU_Table_Program_Chewbacca(ChewbaccaProgram):
    name = "bayeshammer"


    def execute_program(self):
        args = self.args
        self.build_otu_chewbacca(args.outdir, args.groups, args.samples, args.barcodes)



    # TODO should the build table command support parallel operations?
    def build_otu_chewbacca(self, outdir, groups_file, samples_file, barcodes_file):
        """Builds the unannotated OTU table using custom chewbacca script.

        :param outdir: The directory where the matrix should be written.
        :param groups_file: A .groups file containing the OTU names and their consituent/replicant sequences.
        :param samples_file: A .samples file containing the samples that each sequence in the .groups file belongs to.
        :param barcodes_file: A .barcodes file listing all sample names.
        """
        makeDirOrdie(outdir)
        groups = getInputFiles(groups_file)
        debugPrintInputInfo(groups, "read.")

        samples = getInputFiles(samples_file)
        debugPrintInputInfo(samples, "read.")

        barcodes = getInputFiles(barcodes_file)
        debugPrintInputInfo(barcodes, "read.")

        printVerbose("Building matrix...")
        buildOTUtable(groups, samples, barcodes[0], "%s/%s.txt" % (outdir, "matrix"))
        printVerbose("Done building.")