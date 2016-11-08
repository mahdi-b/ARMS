from Bio import SeqIO


class BufferedSeqWriter:
    """An object that handles the buffered writing of seqences to files.  Be sure to call flush() at some point."""

    def __init__(self, output_f, filetype, max_buff=5000):
        self.len = 0
        self.buffer = []
        self.out_stream = open(output_f,'w')
        self.filetype = filetype
        self.max_buff = max_buff

    def write(self, seq):
        self.buffer.append(seq)
        self.len +=1
        if self.len > self.max_buff:
            self.__write__()


    def __write__(self):
        SeqIO.write(self.buffer, self.out_stream, self.filetype)
        self.buffer = []
        self.len = 0


    def flush(self):
        self.__write__()
        self.out_stream.close()