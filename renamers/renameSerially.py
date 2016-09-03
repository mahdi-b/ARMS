from Bio import SeqIO
import sys
import os
from classes.Helpers import getFileName, strip_ixes, clip_count


def serialRename(input_file, output_file, file_type, clip=True):
    """Takes in a fasta file and outputs a new fasta with the sequences renamed.  Renaming convention is x.y.z<n> for
        x.y.z.fasta, where n is an integer in the range [0:n] where n is the position of the sequence in the input_file.
        e.g. The sequences in SiteX_SampleA.fasta are renamed:
                SiteX_SampleA_0, SiteX_SampleA_1, SiteX_SampleA_2, etc.
    :param input_file:      Input fasta or fastq file.
    :param output_file:     Filepath for the output file.
    :param file_type:       "fasta" or "fastq"
    """

    names_file = "%s/%s.names" % (os.path.dirname(output_file), getFileName(input_file))
    seqPrefix = strip_ixes(input_file)
    i = 0
    with open(output_file, 'w') as output:
        with open(names_file,'w') as log:
            for s in SeqIO.parse(input_file, file_type):
                s.id ="%s_ID%s" % (seqPrefix, i)
                log.write("%s\t%s\n" % (s.id, clip_count(seqPrefix, '_')))
                s.description = ""
                SeqIO.write(s, output, file_type)
                i += 1
        log.close()
    output.close()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Usage: input_file output_file renaming_prefix file_type"
    else:
        serialRename(*sys.argv[1:4])