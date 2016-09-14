import pandas as pd

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
        rslt_row_count = int(pct*num_rows)
    else:
        rslt_row_count = n

    return dataframe[:rslt_row_count]


def makeBiome(dataframe):
    """Converts a dataframe into a barebones .biom file.

    :param dataframe: The dataframe to convert
    :return: The text content of the converted dataframe.


    Example .biome file:

    {
        "id": null,
        "format": "Biological Observation Matrix 0.9.1-dev",
        "format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
        "type": "OTU table",
        "generated_by": "QIIME revision 1.4.0-dev",
        "date": "2011-12-19T19:00:00",
        "rows": [
            {"id": "GG_OTU_1", "metadata": null},
            {"id": "GG_OTU_2", "metadata": null},
            {"id": "GG_OTU_3", "metadata": null},
            {"id": "GG_OTU_4", "metadata": null},
            {"id": "GG_OTU_5", "metadata": null}
        ],
        "columns": [
            {"id": "Sample1", "metadata": null},
            {"id": "Sample2", "metadata": null},
            {"id": "Sample3", "metadata": null},
            {"id": "Sample4", "metadata": null},
            {"id": "Sample5", "metadata": null},
            {"id": "Sample6", "metadata": null}
        ],
        "matrix_type": "dense",
        "matrix_element_type": "int",
        "shape": [5, 6],
        "data": [[0, 0, 1, 0, 0, 0],
                 [5, 1, 0, 2, 3, 1],
                 [0, 0, 1, 4, 2, 0],
                 [2, 1, 1, 0, 0, 1],
                 [0, 1, 1, 0, 0, 0]]
    }
    """
    header = "\t\"id\": \"gregs data\",\n\
        \"format\": \"Biological Observation Matrix 0.9.1-dev\",\n\
        \"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\",\n\
        \"type\": \"OTU table\",\n\
        \"generated_by\": \"QIIME revision 1.4.0-dev\",\n\
        \"date\": \"2011-12-19T19:00:00\","
    row_content = ["\t\t{\"id\": \"%s\", \"metadata\": \" \"},"%row for row in dataframe.index.values]
    row_content[-1] = row_content[-1][:-1]
    rows = "\t\"rows\": [\n%s\n\t]," % "\n".join(row_content)

    col_content = ["\t\t{\"id\": \"%s\", \"metadata\": \" \"},"% col for col in dataframe.columns.values]
    col_content[-1] = col_content[-1][:-1]
    cols = "\t\"columns\": [\n%s\n\t]," % "\n".join(col_content)

    predata = "\t\"matrix_type\": \"dense\",\n\
        \"matrix_element_type\": \"int\",\n\
        \"shape\": [%d, %d]," % (len(dataframe.index.values), len(dataframe.columns.values))
    data_content = str(dataframe.values.tolist()).replace('], ','],\n\t\t')
    data = "\t\"data\": %s" % data_content

    body = "%s\n%s\n%s\n%s\n%s" % (header, rows, cols, predata, data)
    text = "{\n%s\n}"% body
    print text
    return text
