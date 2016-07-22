import signal
import os
from Helpers import printVerbose
import re 
import Validators  # fileExists, installed, etc..
import subprocess

# TODO: move programs paths to a config file
MACSE = "/ARMS/programs/MACSE/macse_v1.01b.jar"
FASTX = "~/programs/fastx/bin/"
PEAR = "~/programs/pear-0.9.4-bin-64/pear-0.9.4-64"


class ProgramRunner(object):
    commands = {
        "barcode.splitter": FASTX + "/fastx_barcode_splitter.pl  --bcfile %s --prefix %s --suffix .fastq --bol \
                                --mismatches 1",
        "fastx_renamer": FASTX + "/fastx_renamer -n COUNT -i %s -o %s -Q 33",
        "pear": PEAR + " -f %s -r %s -o %s -j %s",
        "make.contigs": "mothur #make.contigs(ffastq=%s, rfastq=%s, bdiffs=1, pdiffs=2, oligos=%s, processors=%s)",
        "fastq.info": "mothur #fastq.info(fastq=%s)",
        "trim.seqs": "mothur #trim.seqs(fasta=%s, oligos=%s, maxambig=1, maxhomop=8, minlength=300, maxlength=550, \
                            bdiffs=1, pdiffs=2)",
        "align.seqs": "mothur #align.seqs(candidate=%s, template=%s, flip=t)",
        "unique.seqs": "mothur #unique.seqs(fasta=%s)",
        "macse_align": "java -jar " + MACSE + " -prog enrichAlignment  -seq %s -align %s -seq_lr %s -maxFS_inSeq 0 \
                            -maxSTOP_inSeq 0  -maxINS_inSeq 0 -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 \
                            -out_NT %s_NT -out_AA %s_AA -seqToAdd_logFile %s_log.csv",
        "macse_format": "java -jar " + MACSE + "  -prog exportAlignment -align %s -charForRemainingFS - -gc_def 5 \
                        -out_AA %s -out_NT %s -statFile %s",
        "chmimera.uchime": "mothur #chimera.uchime(fasta=%s, name=%s)",
        "remove.seqs": "mothur #remove.seqs(accnos=%s, fasta=%s)",
        "cluster-swarm": "",
        }

    def __init__(self, program, params, conditions={}):
        """Initalizes a ProgramRunner object.

        :param program: A string from the dictionary of available commands.
        :param params: A fully initalized list of parameters for the program to run.
        :param conditions:
        """

        self.program = program
        self.command = self.commands[program] % tuple(params)
        self.conditions = conditions

    def validateConditions(self, conditions):
        """Validates all conditions using a Validator object, returning True if successful and raising an exception
            if not.

        :param conditions: A list of pairs where the first element in each pair is a <Validator.function> and the
                            second element is a list of parameters to that function.
                            i.e. [(<Validator.function>, [parameters to <Validator.function>],) , ]
        :return: True if validation is successful, and raising an exception in <Validator.function> if not.
        """
        # if any of the conditions fail, an exception is raised in getattr
        # TODO: The exception should be raised in this function, so users know that the error is in validation of their conditions, not the
        #   accessing of properties.
        for condition in conditions.iteritems():
            getattr(Validators, condition[0])(condition[1])
        return True


    def dryValidateConditions(self, conditions):
        """Prints validation procedures without actually executing them.

        :param conditions: A list of pairs where the first element in each pair is a <Validator.function> and the
                            second element is a list of parameters to that function.
                            i.e. [(<Validator.function>, [parameters to <Validator.function>],) , ]
        """
        # if any of the conditions fails, stop execution
        # TODO ask Mahdi if this is the intended behavior
        for condition in conditions.iteritems():
            print "\t\tvalidating that %s, %s" % (str(condition[1]), condition[0])
            # TODO: add code

    def splitCommand(self, command):
        """Splits a string on whitespace.  If it is a mothur command, split only on the first space.
        :param command: A fully paramaterized command string
        :return: A command string split on white space (or only on the first white space if it is a mothur command.)
        """
        # expand the ~ in the command line
        command = os.path.expanduser(command)
        commandList = []
        if "mothur" in command:
            commandList = ["mothur",  re.sub(r'mothur\s+', "", command)]
        else:
            commandList = command.split()
        return commandList

    def run(self):
        """Validates conditions (or prints them for a dry run), then executes the command.

        :return: None
        """
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
        """Prints the validation commands that would be performed in an actual run, and the command that would be run.

        :return:
        """
        self.dryValidateConditions(self.conditions)
        return self.command
