import matplotlib
matplotlib.use('Agg')
from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import printVerbose, makeDirOrdie, getInputFiles, strip_ixes
from matplotlib import pyplot as plt
from Visualize_Helpers import subset_dataframe


class Visualize_OTU_Sample_Composition_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"

    def execute_program(self):
        args = self.args
        input_file = getInputFiles(args.input_f)[0]
        output_file = "%s/%s.png" % (args.outdir, strip_ixes(input_file))
        data_frame = subset_dataframe(args.input_f, args)
        self.visualize_otu_sample_comp(data_frame, output_file)

    def visualize_otu_sample_comp(self, data_frame, output_file):
        """Creates a stacked barchart showing the OTU composition in each sample.

        :param data_frame: A pandas dataframe to graph.
        :param output_file: Filepath to the output graphics file.
        """
        ncols = len(data_frame.columns.values)
        nrows = len(data_frame.index.values)
        printVerbose("Computing dataframe values...")
        sums = data_frame.sum(0).values
        data_frame = data_frame.divide(sums)
        printVerbose("Transposing dataframe...")
        data_frame.transpose().plot(kind='bar', stacked=True, ylim=(0, 1), figsize=(10 + ncols / .5, 7 + nrows / 10.0),
                               colormap=plt.cm.hsv)
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=5)
        printVerbose("Saving image %s..." % output_file)
        plt.savefig(output_file, dpi=200)
