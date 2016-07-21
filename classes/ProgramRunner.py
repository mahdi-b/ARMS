import signal
import os
from Helpers import printVerbose
import re 
import Validators as val # fileExists, installed, etc..
import subprocess

# TODO: move programs paths to a config file
MACSE= "/ARMS/programs/MACSE/macse_v1.01b.jar"
FASTX = "~/programs/fastx/bin/"
PEAR= "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64"


class ProgramRunner(object):
    commands = {
        "barcode.splitter": FASTX+"/fastx_barcode_splitter.pl  --bcfile %s --prefix %s --suffix .fastq --bol --mismatches 1",
        "fastx_renamer": FASTX+"/fastx_renamer -n COUNT -i %s -o %s -Q 33",
        "pear": PEAR+" -f %s -r %s -o %s -j %s",
        "make.contigs": "mothur #make.contigs(ffastq=%s, rfastq=%s, bdiffs=1, pdiffs=2, oligos=%s, processors=%s)",
        "fastq.info": "mothur #fastq.info(fastq=%s)",
        "trim.seqs": "mothur #trim.seqs(fasta=%s, oligos=%s, maxambig=1, maxhomop=8, minlength=300, maxlength=550, bdiffs=1, pdiffs=2)",
        "align.seqs": "mothur #align.seqs(candidate=%s, template=%s, flip=t)",
        "unique.seqs": "mothur #unique.seqs(fasta=%s)",
        "macse_align": "java -jar "+MACSE+" -prog enrichAlignment  -seq %s -align %s -seq_lr %s -maxFS_inSeq 0 -maxSTOP_inSeq 0  -maxINS_inSeq 0 -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT %s_NT -out_AA %s_AA -seqToAdd_logFile %s_log.csv",
        "macse_format": "java -jar " +MACSE+"  -prog exportAlignment -align %s -charForRemainingFS - -gc_def 5 -out_AA %s -out_NT %s -statFile %s",
        "chmimera.uchime": "mothur #chimera.uchime(fasta=%s, name=%s)",
        "remove.seqs": "mothur #remove.seqs(accnos=%s, fasta=%s)",
        "cluster-swarm": "",
        }

    def __init__(self, program, params, conditions={}):
        self.program = program
        self.command = self.commands[program] % tuple(params)
        self.conditions = conditions

    def validateConditions(self, conditions):
        # if any of the conditions fails, eception is raised in getattr
        for condition in conditions.iteritems():
            getattr(val, condition[0])(condition[1])
        return True

    def dryValidateConditions(self, conditions):
        # if any of the conditions fails, stop execution
        for condition in conditions.iteritems():
            print "\t\tvalidating that %s, %s" % (str(condition[1]), condition[0])

    def splitCommand(self, command):
        # expand the ~ in the command line
        command = os.path.expanduser(command)
        commandList=[]
        if "mothur" in command:
            commandList = ["mothur",  re.sub(r'mothur\s+', "", command)]
        else:
            commandList = command.split()
        return commandList

    def run(self):
        try:
            self.validateConditions(self.conditions)
            if printVerbose.VERBOSE:
                self.dryRun()
            commandList = self.splitCommand(self.command)
            print commandList
            with open(os.devnull, 'w') as fnull:
                subprocess.call(commandList, stderr=fnull, stdout=fnull)
        except KeyboardInterrupt:
            return 

    def dryRun(self):
        self.dryValidateConditions(self.conditions)
        return self.command
