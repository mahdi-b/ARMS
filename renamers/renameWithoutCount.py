import sys
from Bio import SeqIO
from classes.Helpers import clip_count


def removeCountsFromFastFile(input_file, output_file, file_type, clip_char="_"):
    """Removes the counts from a fasta/fastq file. e.g. "BALI_4462_0_13" becomes "BALI_4462_0".

    :param input_file: Input fasta/fastq file
    :param output_file: Output fasta/fastq file
    :param file_type: Either 'fasta' or 'fastq'
    :param clip_char: The character that delimits the count.  Identifies the suffix to cut at.
                        Defaults to '-'.
    :return: The output filename
    """
    mySeqs=[]
    i = 0
    output_handle = open(output_file, "w")
    for mySeq in SeqIO.parse(input_file, file_type):
        name = clip_count(mySeq.id, clip_char)
        mySeq.id = name
        mySeq.name = ""
        mySeq.description = ""
        mySeqs.append(mySeq)
        if i % 5000 == 0:
            SeqIO.write(mySeqs, output_handle, file_type)
            mySeqs = []
        i += 1
    SeqIO.write(mySeqs, output_handle, file_type)
    return output_file


def removeCountsFromNamesFile(input_names_file, output_names_file, clip_char="_"):
    """Removes the counts all names in a names file.  e.g.
    "BALI_4462_0_ID_13_13  BALI_4462_0_ID_14_5 BALI_4462_0_ID_15_8"
    becomes
    "BALI_4462_0_ID_13  BALI_4462_0_ID_14 BALI_4462_0_ID_15".


    :param input_file: Input names file
    :param output_file: Output names file
    :param clip_char: The character that delimits the count.  Identifies the suffix to cut at.
                        Defaults to '-'.
    :return: The output filename
    """
    with open(output_names_file, "w") as output:
        for line in open(input_names_file, 'r'):
            data = line.split(" ")
            seed_name = clip_char.join(data[0].split(clip_char)[:-1])
            children = []
            for item in data[1:]:
                children.append(clip_char.join(item.split(clip_char)[:-1]))
            output.write("%s\t%s\n" % (seed_name, " ".join(children)))
    return output_names_file