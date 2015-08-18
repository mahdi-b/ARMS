import argparse
import logging
import sys
from multiprocessing import Pool
from classes.Helpers import printVerbose
from classes.Helpers import makeDirOrdie
import os


from classes.ProgramRunner import ProgramRunner

version = "0.01"
FORMAT = "%(asctime)-15s  %(message)s"
args={}

def runInstance(myInstance):
   # Use the global version to facilitate calling workers
   if args.dryRun:
      logging.info(myInstance.dryRun())   
   else:
      logging.info(myInstance.dryRun())   
      myInstance.run()



def preprocessData(args, pool=Pool(processes=1)):
   # TODO: test run.name is a single word


   makeDirOrdie(args.outDir) 



   printVerbose("Preprocessing the data:")

   # *****************************************************************************************
   printVerbose("\t Renaming sequences")
   # "~/programs/fastx/bin/fastx_renamer -n COUNT -i %s %s"
   rename_outFile_f = os.path.join("outDir/", os.path.basename(args.input_f)+"_renamed")
   rename_outFile_r = os.path.join("outDir/", os.path.basename(args.input_r)+"_renamed")

   pool.map(runInstance, [ProgramRunner("fastx_renamer",[args.input_f, rename_outFile_f], {"exists":[args.input_f]}),
                          ProgramRunner("fastx_renamer",[args.input_r, rename_outFile_r], {"exists":[args.input_r]}),
                          ])
   printVerbose("\tRenamed X sequences")
   # *****************************************************************************************
   # Making the contigs using Pear
   # "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s"
   assembledPrefix = os.path.join("outDir", args.name)
   pool.map(runInstance, [ProgramRunner("pear", 
                                           (rename_outFile_f, rename_outFile_r, assembledPrefix, args.threads), 
                                           {"exists":[rename_outFile_f, rename_outFile_r]}) 
                          ])
   assembledFastqFile = os.path.join("outDir", args.name+".assembled.fastq")
   # add py char to a web-page
   printVerbose("\t %s sequences assembled, %s contigs discarded, %s sequences discarded" % (1,1,1))


   # *****************************************************************************************
   # Converting fastq to fasta file (do with mothur or BioSeqIO to keep prog deps to a minimum)
   pool.map(runInstance, [ProgramRunner("fastq.info", 
                                           [assembledFastqFile], 
                                           {"exists": [assembledFastqFile]}) 
                          ])
   assembledFastaFile = os.path.splitext(assembledFastqFile)[0]+".fasta"
   # TODO: add py char to a web-page
   printVerbose("\t converted fastq to fasta")
   # *****************************************************************************************
   # Trimming and assigning reads to groups
   # trim.seqs(fasta=%, oligos=%s, maxambig=0, maxhomop=8, minlength=300, maxlength=550, bdiffs=1, pdiffs=2)

   pool.map(runInstance, [ProgramRunner("trim.seqs",
                                        [assembledFastaFile, args.barcodes],
                                        {"exists": [assembledFastaFile]})
                          ])
   printVerbose("\t %s sequences were assigned to groups and %s sequences were discareded")
   trimmedFasaFile = os.path.splitext(assembledFastqFile)[0]+".trim.fasta"
   # *****************************************************************************************
   # Aligning against the BIOCODETEMPLATE database
   pool.map(runInstance, [ProgramRunner("align.seqs",
                                        [trimmedFasaFile, args.db],
                                        {"exists": [trimmedFasaFile]})
                          ])
   printVerbose("\t %s sequences were assigned to groups and %s sequences were discareded")



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
    parser_preprocess.add_argument('-d', '--db', required=True, help="Db against which the seqeunces are aligned")

    parser_preprocess.set_defaults(func=preprocessData)
    
    global args
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format=FORMAT, level=logging.DEBUG)
    else:
        logging.basicConfig(format=FORMAT, level=logging.ERROR)

    printVerbose.VERBOSE = args.verbose
    printVerbose("Running with %s threads" % args.threads)
    pool = Pool(processes=args.threads)
    logging.debug("Initial ARGS are:")   
    logging.debug(args)
    args.func(args, pool)   

if __name__ == "__main__":
    main(sys.argv)
   


