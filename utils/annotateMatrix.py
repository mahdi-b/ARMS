from classes.Helpers import clip_count


def buildLine(row_name_str, child_list, header_padding, padding):
    out = "%-*s" % (header_padding, row_name_str)
    for child in child_list:
        out += "%*s" % (padding, child)
    out += "\n"
    return out

# TODO clipping?
def annotateMatrix(matrix_file, annotation_file,  out_file, id_col=0, tax_col=5, padding=15,
                   clip_count_from_annotations=True):
    # CREATE A DICTIONARY MAPPING SEQUENCE ID TO TAXONOMIC NAMES
    id_to_tax = {}
    for line in open(annotation_file, 'r'):
        data = line.split()
        id = data[id_col].rstrip()
        if clip_count_from_annotations:
            id = clip_count(id)
        tax = data[tax_col].rstrip()
        id_to_tax[id] = tax
    identifiable = id_to_tax.keys()
    # COMPUTE THE LONGEST OTU NAME IN THE ANNOTATION FILE (FORMATTING PURPOSES)
    max_row_name_length = 0
    for key in identifiable:
        name_len = len(id_to_tax[key])
        if name_len > max_row_name_length:
            max_row_name_length = name_len

    # Parse through the matrix file to see what the current longest row name is
    for line in open(matrix_file, 'r'):
        row_name_length = len(line.split()[0].rstrip())
        if row_name_length > max_row_name_length:
            max_row_name_length = row_name_length

    # Parse through the matrix file again to replace any identifiable ids and reformat the file
    with open(out_file, 'w') as out:
        for line in open(matrix_file, 'r'):
            data = line.split()
            current_id = data[0].rstrip()
            if  current_id in identifiable:
                tax = id_to_tax[current_id]
                out_line = buildLine(tax, data[1:], max_row_name_length, padding)
                out.write(out_line)
            else:
                out_line = buildLine(data[0], data[1:], max_row_name_length, padding)
                out.write(out_line)
