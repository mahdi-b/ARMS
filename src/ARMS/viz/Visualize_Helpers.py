import pandas as pd
import math
from classes.Helpers import *

def subset_dataframe(input_f, args):
    # Gather input files
    inputs = getInputFiles(input_f)
    input_table = inputs[0]
    debugPrintInputInfo(inputs, "visualize")
    data_frame = otu_table_to_dataframe(input_table)
    data_frame = select_top_otus(data_frame, n=100)

    if args.pct is not None:
        printVerbose("Subsetting dataframe with pct = %f" % args.pct)
        if args.pct <= 0:
            print "Error: 'pct' parameter must be a positive real number in the range (0,1].\n"
            exit()
        return select_top_otus(data_frame, pct=args.pct, n=0)

    elif args.count is not None:
        printVerbose("Subsetting dataframe with count = %d" % args.count)
        if args.count < 1:
            print "Error: 'count' parameter must be a positive integer."
            exit()
        return select_top_otus(data_frame, pct=0, n=args.count)

    elif args.names is not None:
        printVerbose("Subsetting dataframe with namesfile = %s" % args.names)
        names_file = getInputFiles(args.names, ignore_empty_files=False)[0]
        return select_otus_by_name_file(data_frame, names_file)

    else:
        return data_frame


def otu_table_to_dataframe(path):
    """Parses an OTU table into a pandas DataFrame, and returns it.  OTU tables should have the following format:
    OTU SampleName  SampleName  SampleName ...
    OTU_Name    1   1   0
    OTU_Name    2   3   4
    ...
    i.e. the first row should be a header column, with the first word being "OTU", followd by a tab and a tab-delimited
    list of SampleNames/locations.
    Each subsequent row should contain an OTU name, a tab, and a tab-delimited list of the OTU's abundance in/at each
    Sample.

    :param path: Filepath to the white-space delimited OTU matrix.
    :return: A pandas dataframe where dataframe.columns returns the sample names, and dataframe.index returns the
            OTU names.
    """
    frame = pd.read_table(path, delim_whitespace=True, index_col=0, header=0)
    return frame


def select_otus_by_name(dataframe, names):
    """Given a dataframe, and a list of OTU names, returns a dataframe with just the indicated names.

    :param dataframe: The dataframe to pull data from.
    :param names: A list of OTU names to include in the returned dataframe.
    :return: A dataframe with just the named OTUs.
    """
    return dataframe.loc[dataframe.index.isin(names)]


def select_otus_by_name_file(dataframe, names_filepath):
    """Given a dataframe, and a file containing a list of OTU names, returns a dataframe with just the indicated names.
        Each line of the specified file should contain a single OTU name to pull.
        e.g.

        OTU_Name_3
        OTU_Name_2
        OTU_Name_5
        ...

    :param dataframe: The dataframe to pull data from.
    :param names_filepath: Filepath to a file of a list of desired OTU names.
    :return: A dataframe with just the named OTUs.
    """
    names = [line.strip() for line in open(names_filepath, 'r')]
    return select_otus_by_name(dataframe, names)


def select_top_otus(dataframe, pct=.1, n=0):
    """Selects the most abundant OTUs from a dataframe, and returns them in a new dataframe.  The number of OTUs to
    return can be specified as either a decimal percentage in the range [0,1], or a specific number of OTUs to return.
    Note: Uses pct as the quantity specifier if n < 1.  Uses n as the quantity specifier if n >= 1.
    :param dataframe: The dataframe to pull OTUs from.
    :param pct: A decimal percentage in the range [0,1] indicating what percentage of the top scoring OTUs should be
                    returned.  e.g. pct = .1 specifies that the top 10% of OTUs (by abundance) should be returned.
                    Default: 0.1.
    :param n: An integer specifying the exact number of OTUs to return.  e.g n=10 returns the 10 most abundant OTUs.
    :return: A dataframe with the specified number of OTUS, sorted by abundance.
    """
    dataframe['Total'] = dataframe.sum(1)
    dataframe = dataframe.sort_values("Total", ascending=False)
    num_rows = len(dataframe.index.values)
    rslt_row_count = num_rows
    if n < 1:
        rslt_row_count = int(math.ceil(pct*num_rows))
    else:
        rslt_row_count = n

    return dataframe[:rslt_row_count]
