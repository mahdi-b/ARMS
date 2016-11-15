import sys
from BufferedWriter import BufferedFileWriter
class OutputManager(object):
    def __init__(self, stdout=sys.stdout, log="log.txt"):
        self.stdout = stdout
        self.log = log
        self.verbosity = debug
        self.writer= BufferedFileWriter(self.log)


    def debug(self, message):
        pass


    def info(self, message):
        pass


    def error(self, message):
        pass


    def log(self, message):
        log.write(message)
