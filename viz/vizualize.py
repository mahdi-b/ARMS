import pandas as pd
import sys
from matplotlib import pyplot as plt
import seaborn as sns
from utils import otu_table_to_dataframe
import numpy as np


def viz_OTU_composition_barchart(frame, output_file="", as_pct_composition=True):
    ncols = len(frame.columns.values)
    nrows = len(frame.index.values)
    if as_pct_composition:
        sums = frame.sum(0).values
        frame = frame.divide(sums)
    frame.transpose().plot(kind='bar', stacked=True, ylim=(0,1), figsize=(ncols/3.0+6, nrows/10.0), colormap=plt.cm.hsv)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=5)
    if output_file != "":
        plt.savefig(output_file, dpi=200)
    else:
        plt.show()


def viz_OTU_abundance_heatmap(frame, output_file=""):

    ncols = len(frame.columns.values)
    nrows = len(frame.index.values)
    fig, ax = plt.subplots(figsize=( 5+ncols/10.0, nrows/10.0))
    heatmap = ax.pcolor(frame, cmap=plt.cm.gnuplot2)
    ax.set_xticks(np.arange(ncols) + 0.5)
    ax.set_yticks(np.arange(nrows) + 0.5)
    ax.set_xticklabels(frame.columns.values,rotation=90, fontsize=4)
    ax.set_yticklabels(frame.index.values, fontsize=4)
    cbar = plt.colorbar(heatmap)
    if output_file != "":
        plt.savefig(output_file, dpi=200)
    else:
        plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: function   input_otu_table   [output_file]"
    else:

        outfile = ""
        option = sys.argv[1]
        OTU_matrix_file_path = sys.argv[2]
        if len(sys.argv) >= 4:
            outfile = sys.argv[3]

        frame = otu_table_to_dataframe(OTU_matrix_file_path)
        frame = select_top_otus(frame, pct=.1)
        print frame
        """
        function = ""
        if option == "heatmap":
            viz_OTU_abundance_heatmap(otu_table, outfile)
        if option == "barchart":
            viz_OTU_composition_barchart(otu_table, outfile)
        """
#path = "testMatrix.txt"

#path = "/home/greg/ARMS/testARMS/16_annotate_bold/matrix.txt"

#viz_OTU_abundance_heatmap(path)