import os
from Helpers import printVerbose

import Validators as val # fileExists, installed, etc..





class ProgramRunner(object):


    commands = {
        "fastx_renamer": "~/programs/fastx/bin/fastx_renamer -n COUNT -i %s -o %s -Q 33",
        "pear": "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64 -f %s -r %s -o %s -j %s",
        "make.contigs": "mothur \"#make.contigs(ffastq=%s, rfastq=%s, bdiffs=1, pdiffs=2, oligos=%s, processors=%s)\"",
        "fastq.info": "mothur \"#fastq.info(fastq=%s)\"",
        "trim.seqs": "mothur \"#trim.seqs(fasta=%s, oligos=%s, maxambig=1, maxhomop=8, minlength=300, maxlength=550, bdiffs=1, pdiffs=2)\"",
        "align.seqs": "mothur \"#align.seqs(candidate=%s, template=%s, flip=t)\"",
        }

    def __init__(self, program, params, conditions={}):
        self.program = program
        self.command = self.commands[program] % tuple(params)
        self.conditions = conditions

    def validateConditions(self, conditions):
        # if any of the conditions fails, eception is raised in getattr
        print "Conditions: %s" % conditions
        for condition in conditions.iteritems():
            getattr(val, condition[0])(condition[1])
        return True

    def dryValidateConditions(self, conditions):
        # if any of the conditions fails, stop execution
        print "Dry-validating the conditions"
        for condition in conditions.iteritems():
            print "validating that %s, %s" % (str(condition[1]), condition[0])

    def run(self):
        print "in run conditions: %s" % self.conditions
        self.validateConditions(self.conditions)
        if printVerbose.VERBOSE:
            self.dryRun()
        os.system(self.command)

    def dryRun(self):
        self.dryValidateConditions(self.conditions)
        return self.command
