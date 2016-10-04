from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import *
from viz.utils import *
from matplotlib import pyplot as plt

class Visualize_OTU_Sample_Composition_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"

    def execute_program(self):
        args = self.args
        self.visualize_otu_sample_comp(args.input_f, args.outdir)

    def visualisze_otu_sample_comp(self, input_f, outdir):
        """Creates a stacked barchart showing the OTU composition in each sample.

        :param input_f: Filepath to an input file or folder to rename.
        :param outdir: Filepath to the output directory.
        """
        # Make the output directory, or abort if it already exists
        makeDirOrdie(outdir)

        # Gather input files
        inputs = getInputFiles(input_f)
        input_table = inputs[0]
        debugPrintInputInfo(inputs, "visualize")
        data_frame = otu_table_to_dataframe(input_table)
        data_frame = select_top_otus(data_frame, n=100)
        output_file = "%s/%s.png" % (outdir, strip_ixes(input_table))
        ncols = len(data_frame.columns.values)
        nrows = len(data_frame.index.values)
        sums = data_frame.sum(0).values
        data_frame = data_frame.divide(sums)
        data_frame.transpose().plot(kind='bar', stacked=True, ylim=(0, 1), figsize=(10 + ncols / .5, 7 + nrows / 10.0),
                               colormap=plt.cm.hsv)
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=5)
        plt.savefig(output_file, dpi=200)
