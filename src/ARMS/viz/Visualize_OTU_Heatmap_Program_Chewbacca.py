import numpy as np
from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import *
from matplotlib import pyplot as plt
from Visualize_Helpers import *

class Visualize_OTU_Heatmap_Program_Chewbacca(ChewbaccaProgram):
    name = "chewbacca"

    def execute_program(self):
        args = self.args
        makeDirOrdie(args.outdir)
        input_file = getInputFiles(args.input_f)[0]
        output_file = "%s/%s.png" % (args.outdir, strip_ixes(input_file))
        data_frame = subset_dataframe(args.input_f, args)
        self.visualize_otu_heatmap(data_frame, output_file)

    def visualize_otu_heatmap(self, data_frame, output_file):
        """Visualizes an OTU table as a heatmap, showing OTU abundance in each Sample.

        :param data_frame: A pandas dataframe to graph.
        :param output_file: Filepath to the output graphics file.
        """
        ncols = len(data_frame.columns.values)
        nrows = len(data_frame.index.values)
        printVerbose("Computing dataframe values...")
        fig, ax = plt.subplots(figsize=(10 + ncols / .5, 7 + nrows / 10.0))
        heatmap = ax.pcolor(data_frame, cmap=plt.cm.gnuplot2)
        ax.set_xticks(np.arange(ncols) + 0.5)
        ax.set_yticks(np.arange(nrows) + 0.5)
        ax.set_xticklabels(data_frame.columns.values, rotation=90, fontsize=4)
        ax.set_yticklabels(data_frame.index.values, fontsize=4)
        cbar = plt.colorbar(heatmap)
        printVerbose("Saving image %s..." % output_file)
        plt.savefig(output_file, dpi=200)
