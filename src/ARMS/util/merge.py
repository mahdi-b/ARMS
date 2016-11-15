from classes.BufferedWriter import BufferedFileWriter

def merge_files(input_file_list, output_file):
    """Concatenates the contents of input_file_list into one output_file.  Overwrites any pre-existing content in
        output_file.

    :param input_file_list: List of filepaths to join.
    :param output_file: File path to the output file.  Pre-existing contents will be overwritten
    """
    merged_file = BufferedFileWriter(output_file)
    for fname in input_file_list:
        for line in open(fname):
            merged_file.write(line)
    merged_file.flush()
    return output_file