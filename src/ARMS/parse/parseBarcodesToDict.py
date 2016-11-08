
def parse_barcodes_to_dict(barcodes_file):
    """Parses a two-column .barcodes file to a dictionary mapping sample names to barcode sequences.

    :param barcodes_file: The two-column .barcodes file to parse
    :return: A dictionary mapping sample names to barcode sequences.
    """
    with open(barcodes_file, 'r') as input_f:
        rslt = {}
        for line in input_f:
            data = map(str.rstrip, line.split())
            if len(data) == 2:
                rslt[data[0]] = data[1]
        return rslt