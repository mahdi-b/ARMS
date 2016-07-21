import os
from Bio import SeqIO
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
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



def splitFileBySample(fastaFile, groupsFile, splitFastaDir):
   # Reads the sequences/group association from the group file (similar to mothur)
   # and outputs as many files as there are samples
   # not ideal parallelization since some samples might be larger than others
   # however, this is ideal require sample file independently
   # this also sanitizes the sequences by removing "." dots added by align.seq
   
   # TODO: test that the file fits into memory, otherwise this could cause problems
   mySequences = SeqIO.to_dict(SeqIO.parse(fastaFile, 'fasta'))

   seqsByGroup=defaultdict(list)
   for line  in open(groupsFile, 'r'):
      seq, group = line.rstrip().split()
      seqsByGroup[group].append(seq)
   for group in seqsByGroup.keys():
      outFile = open(os.path.join(splitFastaDir, group), 'a')
      for seq in seqsByGroup[group]:
         SeqIO.write(SeqRecord(mySequences[seq].seq.ungap(".").ungap("-"),id=seq, description=""), outFile, 'fasta')
      outFile.close()
   
      
