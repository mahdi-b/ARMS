import glob
import argparse
import logging
import sys
from multiprocessing import Pool
from classes.Helpers import printVerbose
from classes.Helpers import makeDirOrdie
from classes.Helpers import splitFileBySample
import os
import signal
import time 
from Bio import SeqIO
from Bio.Seq import Seq
from shutil import copyfile


 
from classes.ProgramRunner import ProgramRunner

version = "0.01"
FORMAT = "%(asctime)s  %(message)s"
DATEFMT= "%m/%d %H:%M:%S"


args={}

def runInstance(myInstance):
   # Use the global version to facilitate calling workers
   if args.dryRun:
      logging.info(myInstance.dryRun())   
   else:
      logging.info(myInstance.dryRun())   
      myInstance.run()
      print "-------" * 10



def preprocessData(args, pool=Pool(processes=1)):
   # TODO: test run.name is a single word
   makeDirOrdie(args.outDir) 
   printVerbose("Preprocessing the data:")

   # *****************************************************************************************
   printVerbose("\tRenaming sequences")

   # "~/programs/fastx/bin/fastx_renamer -n COUNT -i %s %s"
   rename_outFile_f = os.path.join(args.outDir, os.path.basename(args.input_f)+"_renamed")
   rename_outFile_r = os.path.join(args.outDir, os.path.basename(args.input_r)+"_renamed")

   # pool.map(runInstance, [ProgramRunner("fastx_renamer",[args.input_f, rename_outFile_f], {"exists":[args.input_f]}),
   #                        ProgramRunner("fastx_renamer",[args.input_r, rename_outFile_r], {"exists":[args.input_r]}),
   #                        ])
   printVerbose("\tRenamed %s sequences")
   # *****************************************************************************************
   # Making the contigs using Pear
   # "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s"
   assembledPrefix = os.path.join(args.outDir, args.name)
   # pool.map(runInstance, [ProgramRunner("pear", 
   #                                         [rename_outFile_f, rename_outFile_r, assembledPrefix, args.threads], 
   #                                         {"exists":[rename_outFile_f, rename_outFile_r]}) 
   #                        ])
   assembledFastqFile = os.path.join(args.outDir, args.name+".assembled.fastq")
   printVerbose("\t%s sequences assembled, %s contigs discarded, %s sequences discarded" % (-1,-1,-1))

   # MAHDI  Comment the rest of the program to replace the pipeline with 
   # 1- split_libraries_fastq.py
   # 2- usearch -fastq_filter
   printVerbose("Splitting based on barcodes")
   pool.map(runInstance, [ProgramRunner("barcode.splitter",
                                           [assembledFastqFile, args.barcodes,  os.path.join(args.outDir,"splitOut_")],
                                           {"exists": [assembledFastqFile]})
                          ])

   listOfSamples =  glob.glob(os.path.join(args.outDir,"splitOut_*"))



def splitFile(args, pool):
   # Split the cleaned file resulting from preprocessData into a user defined
   # number of chunks
   makeDirOrdie(args.outDir)
   splitFileBySample(args.inputFasta, args.groups, args.outDir)
   printVerbose("\tDone splitting file")

# TODO: eventually send a param to Program running, prevenint it from starting after CTRL+C has been invoked
def clean(args, pool):
   makeDirOrdie(args.outDir)
   try:
      if args.program == "macse":
         printVerbose("\t %s Aligning reads using MACSE")

         pool.map(runInstance, [ProgramRunner("macse_align", 
                                                    [ args.db, args.db, os.path.join(args.samplesDir,sample)]+[os.path.join(args.outDir, sample)]*3  
                                                    , {"exists": []}) for sample in os.listdir(args.samplesDir) ])


         printVerbose("\t %s Processing MACSE alignments")
         pool.map(runInstance, [ProgramRunner("macse_format", 
                                                    [os.path.join(args.samplesDir,sample),  os.path.join(args.outDir,sample+"_AA_macse.fasta"), 
                                                     os.path.join(args.outDir,sample+"_NT_macse.fasta"), os.path.join(args.outDir,sample+"_macse.csv")]  
                                                    , {"exists": []}) for sample in os.listdir(args.samplesDir) ])


      pool.close()
      pool.join()

      printVerbose("\tCleaning MACSE alignments")

      # Remove the reference sequences from the MACSE files and remove the non nucleotide characters from the sequences.
      # we need the datbase seq. names to remove them from the results files

      # TODO:IMPORTANT: Merge the files before doing this.
      
      dbSeqNames = SeqIO.to_dict(SeqIO.parse(args.db, "fasta")).keys()
      good_seqs=[]
      samplesList = os.listdir(args.samplesDir)
      print "Will be processing %s samples " % len(samplesList)
      i = 0 
      for sample in  samplesList:
         nt_macse_out = os.path.join(args.outDir,sample+"_NT_macse.fasta")
         for mySeq in SeqIO.parse(nt_macse_out, 'fasta'):
            if mySeq.id not in dbSeqNames:
               mySeq.seq = Seq(str(mySeq.seq[2:]).replace("-", "")) # remove the !! from the beginning
               good_seqs.append(mySeq)
         print "completed %s samples" % i
         i += 1
      SeqIO.write(good_seqs, open( os.path.join(args.outDir,"MACSEOUT_MERGED.fasta"), 'w'), 'fasta')

      printVerbose("\t%s sequences cleaned, %s sequences retained, %s sequences discarded" % (1, 1,1))
   except KeyboardInterrupt:
        pool.terminate()


def removeChimeras(args, pool):
   if args.program == "uchime":

      # *****************************************************************************************
      # findind the chimeras
      # "chmimera.uchime": "mothur #chimera.uchime(fasta=%s, name=%s)",
      pool.map(runInstance, [ProgramRunner("chmimera.uchime", 
                                           [args.inputFile, args.namesFile], 
                                           {"exists": [args.inputFile, args.namesFile]})
                             ])
      # *****************************************************************************************
      # removing the chimeras from input file
      # remove.seqs(accnos=Nucleotidealignment.uchime.accnos, fasta=Nucleotidealignment.fasta)

      uchimeAccnos = glob.glob(os.path.join(os.path.dirname(args.inputFile),"*uchime.chimeras"))[0]
      print "outFile is %s " % uchimeAccnos
      p = ProgramRunner("remove.seqs", [ uchimeAccnos, args.inputFile], {"exists": [uchimeAccnos]})
      runInstance(p)

      # *****************************************************************************************
      # Renaming the outfile and updating the names file
      # remove from the names file the sequences that were removed
      
      pickOutFile = glob.glob(os.path.join(os.path.dirname(args.inputFile),"*pick.fasta"))[0]
      copyfile(pickOutFile, args.outFasta)

      # TO CONITNUE: UPDATE THE NAMES file 
      

   else:
      raise Exception("unknown program %s for chimera detection or removal" % args.program)


   printVerbose("\t removed %s chimeric sequences")


def dropShort(args, pool):
   # Dropping short sequences
   good_seqs =[]
   for seq in SeqIO.parse(args.inputFasta, "fasta"):
      if len(seq.seq) >= int(args.minLenght):
         good_seqs.append(seq)
      else:
         print "seq %s too short (%s bases)" % (seq.id, len(seq.seq))


      



   

def main(argv):
    parser = argparse.ArgumentParser(description="arms description", epilog="arms long description")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+version)
    parser.add_argument("--verbose", default=True, help="increase output verbosity", action="store_true")   
    parser.add_argument('-t', '--threads', type=int, default = 1)
    parser.add_argument('--dryRun', action='store_true', default = False)
    subparsers = parser.add_subparsers(dest='action', help='Available commands')

    # preprocess data
    parser_preprocess = subparsers.add_parser('preprocess')
    parser_preprocess.add_argument('-n', '--name', required=True, help="Run Id")
    parser_preprocess.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads")
    parser_preprocess.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads")
    parser_preprocess.add_argument('-b', '--barcodes', required=True, help="Tab delimted files of barcodes and their samples")
    parser_preprocess.add_argument('-o', '--outDir', required=True, help="Directory where outputs will be saved")
    parser_preprocess.set_defaults(func=preprocessData)
    
    # split file by samples
    parser_split = subparsers.add_parser('partition')
    parser_split.add_argument('-n', '--name', required=True, help="Run Id")
    parser_split.add_argument('-i', '--inputFasta', required=True, help="Input fasta file to split")
    parser_split.add_argument('-g', '--groups', required=True, help=" Groups file in same format as that generated by mothur")
    parser_split.add_argument('-o', '--outDir', required=True, help="Output directory")
    parser_split.set_defaults(func=splitFile)    


    # Clean files
    parser_clean = subparsers.add_parser('clean')
    parser_clean.add_argument('-n', '--name', required=True, help="Run Id")
    parser_clean.add_argument('-s', '--samplesDir', required=True, help="Directory containing the samples file required for clustering")
    parser_clean.add_argument('-p', '--program', required=True, help="Program name for clenaing. Available options are: ... ")
    parser_clean.add_argument('-d', '--db', required=True, help="Database against which to align and filter reads")
    parser_clean.add_argument('-f', '--namesFile', required=False, help=" Mothus .names file to update", default=None)
    parser_clean.add_argument('-o', '--outDir', required=True, help=" Output directory")
    parser_clean.set_defaults(func=clean)    


    # Remove Chimeras
    parser_chimera = subparsers.add_parser('removeChimeras')
    parser_chimera.add_argument('-n', '--name', required=True, help="Run Id")
    parser_chimera.add_argument('-i', '--inputFile', required=True, help="Clean inputs File")
    parser_chimera.add_argument('-p', '--program', required=False, default="uchime", help="Program for detecting and removing chimeras. Default is uchime")
    parser_chimera.add_argument('-f', '--namesFile', required=True, help="Updated names file")
    parser_chimera.add_argument('-o', '--outFasta', required=True, help="Updated names file")
    parser_chimera.add_argument('-s', '--outNames', required=True, help="Updated names file")
    parser_chimera.set_defaults(func=removeChimeras)    


    # Drop short reads
    parser_dropShort = subparsers.add_parser('dropShort')
    parser_dropShort.add_argument('-n', '--name', required=True, help="Run Id")
    parser_dropShort.add_argument('-i', '--inputFasta', required=True, help="Clean inputs File")
    parser_dropShort.add_argument('-f', '--namesFile', required=True, help="Updated names file")
    parser_dropShort.add_argument('-l', '--minLenght', required=True, help="Min. length to keep")
    parser_dropShort.add_argument('-o', '--outFasta', required=True, help="Output file filtered on lenght")
    parser_dropShort.add_argument('-s', '--outNames', required=True, help="Updated names file")
    parser_dropShort.set_defaults(func=dropShort)    


    # Drop short reads
    parser_cluster = subparsers.add_parser('cluster-swarm')
    #parser_cluster.set_defaults(func=clusterReads)    




    global args
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format=FORMAT, level=logging.DEBUG, datefmt=DATEFMT)
    else:
        logging.basicConfig(format=FORMAT, level=logging.ERROR, datefmt=DATEFMT)



    printVerbose.VERBOSE = args.verbose
    printVerbose("Running with %s threads" % args.threads)
    pool = Pool(args.threads)
    logging.debug("Initial ARGS are: %s", args)   
    print("\t\t")   
    
    args.func(args, pool)   

if __name__ == "__main__":
    main(sys.argv)
   


