from Bio import SeqIO


class BufferedFileWriter(object):
    """An object that handles the buffered writing of text to files.  Be sure to call flush() at some point."""

    def __init__(self, output_f, filetype='fasta', max_buff=5000):
        try:
            self.out_stream = open(output_f,'w')
        except:
            print "ERROR: Could not open the file %s for writing." % output_f
        self.len = 0
        self.buffer = []
        self.filetype = filetype
        self.max_buff = max_buff


    def write(self, data):
        self.buffer.append(data)
        self.len +=1
        if self.len > self.max_buff:
            self.__write__()

    def __write__(self):
        self.out_stream.write("\n".join(self.buffer))
        self.buffer = []
        self.len = 0


    def flush(self):
        self.__write__()
        self.out_stream.close()


class BufferedSeqWriter(BufferedFileWriter):
    """An object that handles the buffered writing of seqences to files.  Be sure to call flush() at some point."""

    def __write__(self):
        SeqIO.write(self.buffer, self.out_stream, self.filetype)
        self.buffer = []
        self.len = 0


