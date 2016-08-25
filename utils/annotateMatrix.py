
def buildLine(row_name_str, child_list):
    out = "%*s" % (header_padding, row_name_str)
    for child in child_list:
        out += "%*s" % (padding, child)
    return out

# TODO clipping?
def annotateMatrix(matrix_file, annotation_file, id_col, tax_col, out_file, clip=True):
    # CREATE A DICTIONARY MAPPING SEQUENCE ID TO TAXONOMIC NAMES
    id_to_tax = {}
    for line in open(annotation_file, 'r'):
        data = line.split("\t")
        id = data[id_col].rstrip()
        tax = data[tax_col].rstrip()
        id_to_tax[id] = tax
    identifiable = id_to_tax.keys()

    # COMPUTE THE LONGEST OTU NAME IN THE ANNOTATION FILE (FORMATTING PURPOSES)
    max_row_name_length = 0
    for key in identifiable:
        key_len = len(key)
        if key_len > max_row_name_length:
            max_row_name_length = key_len

    # Parse through the matrix file to see what the current longest row name is
    for line in open(matrix_file, 'r'):
        row_name_length = len(line.split()[0].rstrip())
        if row_name_length > max_row_name_length:
            max_row_name_length = row_name_length

    # Parse through the matrix file again to replace any identifiable ids and reformat the file
    with open(out_file, 'r') as out:
        for line in open(matrix_file, 'r'):
            data = line.split()
            current_id = data[0].rstrip()
            if  current_id in identifiable:
                tax = id_to_tax[current_id]
                out_line = buildLine(tax, data[1:])
                out.write(out_line)
            else:
                out_line = buildLine(data[0], data[1:])
                out.write(out_line)
