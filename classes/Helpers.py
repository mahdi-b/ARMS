import os



class printVerbose(object):
   def __init__(self, msg, newline =True):
      if printVerbose.VERBOSE:
         if newline:
            print msg
         else:
            print msg,

def makeDirOrdie(dirPath):
   if not os.path.isdir(dirPath):
      os.makedirs(dirPath)
   else:
      pass
      #TODO; Uncomment after testing done
      #logging.error("Split fasta directory %s already exists " % hmmerOutputDir) 
      #sys.exit()
   return dirPath


