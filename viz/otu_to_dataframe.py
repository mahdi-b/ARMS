import pandas as pd

def otu_table_to_dataframe(path):
    frame = pd.read_table(path, delim_whitespace=True, index_col=0, header=0)
    return frame
