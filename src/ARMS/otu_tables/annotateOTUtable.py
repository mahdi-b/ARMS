from itertools import product
from parsers.parseVSearchoutForTaxa import *
from classes.ProgramRunner import *
from classes.Helpers import *




def buildLine(row_name_str, child_list):
    """Given an row name and a list of children, writes the row name, a tab, and then the space-delimited child list.

    :param row_name_str:
    :param child_list:
    :return: A string containing the row name, a tab, and then the space-delimited child list.
    """
    out = "%s" % row_name_str
    for child in child_list:
        out += "\t%s" % child
    out += "\n"
    return out


def annotateOTUtable(otu_file, annotation_file, out_file, id_col=0, tax_col=5, clip_count_from_annotations=True):
    """Given the best hits from a cleaned, annotated (with taxonomic names), Vsearch out file, renames each
        sequence ID in OTU table with its taxonomic name in the Vsearch outfile.

    :param otu_file: The otu file to annotate.
    :param annotation_file: Vsearch output file, cleaned by parsers.ParseVsearchOutForTaxa.py
    :param out_file: Filepath to write the resulting annotated OTU file.
    :param id_col: The column number (zero-indexed) in the Vsearch file containing sequence fasta IDs.
    :param tax_col: The column number (zero-indexed) in the Vsearch file containing the taxonomic IDs.
    :param clip_count_from_annotations: If True, clip dereplication counts from Vsearch sequence fasta IDs.
        Default: True.
    :return: Filepath to the resulting annotated OTU table.
    """
    # CREATE A DICTIONARY MAPPING SEQUENCE ID TO TAXONOMIC NAMES
    id_to_tax = {}
    for line in open(annotation_file, 'r'):
        data = line.split()
        id = data[id_col].rstrip()
        if clip_count_from_annotations:
            id = clip_count(id,'_')
        tax = data[tax_col].rstrip()
        id_to_tax[id] = tax
    identifiable = id_to_tax.keys()

    # Parse through the matrix file again to replace any identifiable ids and reformat the file
    with open(out_file, 'w') as out:
        for line in open(otu_file, 'r'):
            data = line.split()
            current_id = data[0].rstrip()
            if  current_id in identifiable:
                tax = id_to_tax[current_id]
                out_line = buildLine(tax, data[1:])
                out.write(out_line)
            else:
                out.write(line)
    return out_file