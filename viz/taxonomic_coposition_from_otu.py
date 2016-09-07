import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from otu_to_dataframe import otu_table_to_dataframe
import numpy as np


def viz_OTU_composition_barchart(OTU_matrix_file_path, output_file="", as_pct_composition=True):
    frame = otu_table_to_dataframe(OTU_matrix_file_path)
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


def viz_OTU_abundance_heatmap(OTU_matrix_file_path, output_file=""):
    frame = otu_table_to_dataframe(OTU_matrix_file_path)
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


#path = "testMatrix.txt"
path = "mat.txt"
#path = "/home/greg/ARMS/testARMS/16_annotate_bold/matrix.txt"
viz_OTU_abundance_heatmap(path)