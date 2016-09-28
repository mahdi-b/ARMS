

def joinFiles(input_file_list, output_file):
    """Concatenates the contents of input_file_list into one output_file.  Overwrites any pre-existing content in
        output_file.

    :param input_file_list: List of filepaths to join.
    :param output_file: File path to the output file.  Pre-existing contents will be overwritten
    """
    with open(output_file, 'w') as outfile:
        i = 0
        output = ""
        for fname in input_file_list:
            with open(fname) as infile:
                for line in infile:
                    i += 1
                    output += line
                    if i%50000 == 0:
                        outfile.write(output)
                        output = ""
        outfile.write(output)
    outfile.close()
    return output_file