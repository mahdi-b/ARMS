import ConfigParser
import logging
import os
import re
import subprocess
import Validator  # fileExists, installed, etc..
from Helpers import printVerbose, validate
from pipes import quote


class ProgramRunner(object):
    """A class to interact with external command line programs.  The class contains a dictionary of formatted command
    strings.  The class supports validation, sanitization, and debugging of user-supplied parameters.

    Attributes:
        DEFAULT_CONFIG_FILEPATH     The filepath to a config file containing user-specified overrides (such as file
                                        paths to exectuables.
        configsLoaded               Specifies that all config variables have been loaded, and that the file should not
                                        be read again durring the execution of this instance.
        programPaths                Specifies the default location of executibles used in the commands dictionary.
        commandTemplates            A dictionary mapping chewbacca commands to un-paramaterized command line strings.
                                        Used as templates for commands.
    """
    DEFAULT_CONFIG_FILEPATH = "chewbacca.cfg"
    configsLoaded = False
    commandTemplates = {}
    programPaths = {
        "MACSE": "/ARMS/programs/MACSE/macse_v1.01b.jar",
        "FASTX": "~/programs/fastx/bin/",
        "PEAR": ""
    }

    def __init__(self, program, params, conditions={}, stdin=open(os.devnull, 'r'), stdout=open(os.devnull, 'w'),
                 stderr=open(os.devnull, 'w')):
        """Initalizes a ProgramRunner object.  Reads chewbacca.cfg and loads configuration settings, builds the command
        string, and configures stdIO pipes.

        :param program:     A chewbacca command (as a string).  Used as a key to fetch the command string from
                                ProgramRunner.commandTemplates.
        :param params:      A list of parameters for the chosen chewbacca command.
        :param conditions:  A dictionary mapping <Validator>.function names to a list of parameters.  The dictionary
                                key/function name is called on each parameter in the corresponding list of parameters.
                                All parameters must satisfy their <Validator>.function in order for the specified
                                program to execute.
        :param stdin:       A handle to the stdin pipe for the command.  Defaults to os.devnull.
        :param stdout:      A handle to the stdout pipe for the command.  Defaults to os.devnull.
        :param stderr:      A hande to the stderr pipe for the command.  Defaults to os.devnull.
        """
        if not ProgramRunner.configsLoaded:
            ProgramRunner.loadConfigs(self)
            ProgramRunner.initalizeCommands(self)
            ProgramRunner.configsLoaded = True
        logging.debug(params)
        self.program = program
        self.command = self.commandTemplates[program] % tuple(params)
        self.conditions = conditions
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr

    def validateConditions(self, conditions):
        """Validates this program's conditions for execution.

        :param conditions: A dictionary mapping <Validator>.function names to a list of parameters.  The dictionary
                            key/function name is called on each parameter in the corresponding list of parameters.
                            All parameters must satisfy their <Validator>.function in order for the specified
                            program to execute.

        :return:            True if validation is successful, and raising an exception in <Validator.function> if not.
        """
        return validate(conditions)


    def dryValidateConditions(self, conditions):
        """Prints validation procedures without actually executing them.

        :param conditions: A dictionary mapping <Validator>.function names to a list of parameters.  The dictionary
                                key/function name is called on each parameter in the corresponding list of parameters.
                                All parameters must satisfy their <Validator>.function in order for the specified
                                program to execute.
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
            commandList = ["mothur", re.sub(r'mothur\s+', "", command)]
        else:
            commandList = command.split()
        return commandList

    def run(self):
        """Validates conditions (or prints them for a dry run), sanitizes parameters, and then executes the command.

        :return: None
        """
        try:
            self.validateConditions(self.conditions)
            if printVerbose.VERBOSE:
                self.dryRun()
            # call and check_call are blocking, Popen is non-blocking
            print "running " + self.command
            subprocess.Popen(self.command, shell=True)
            # commandList = self.splitCommand(self.command)
            # print commandList
            # print " ".join(commandList)
            # subprocess.call(commandList, stderr=fnull, stdout=fnull)
            # subprocess.call(commandList,stdin=self.stdin,stdout=self.stdout,stderr=self.stderr)
            # subprocess.call(" ".join(commandList))
        except KeyboardInterrupt:
            return

    def dryRun(self):
        """Prints the validation procedures that would be performed, and the commands that would be run in an actual
            run, without actually executing them.

        :return:
        """
        self.dryValidateConditions(self.conditions)
        return self.command

    def loadConfigs(self):
        """Loads configuration settings from the chewbacca configuration file (located in
            ProgramRunner.DEFAULT_CONFIG_FILEPATH).  If the file is not found, default values are used. Overwrites the
            default values of any found entries.  Currently loads the [Program  Paths] section, updating the
            programPaths dictionary.
        """
        configFilePath = ProgramRunner.DEFAULT_CONFIG_FILEPATH
        configSection = "Program Paths"
        if os.path.isfile(configFilePath):
            logging.debug("Loaded Chewbacca config files from %s" % configFilePath)
            config = ConfigParser.RawConfigParser()
            config.read(configFilePath)

            for program in ProgramRunner.programPaths.keys():
                if config.has_option(configSection, program):
                    configSetting = config.get(configSection, program)
                    ProgramRunner.programPaths[program] = configSetting
                    logging.debug("Read %s filepath as %s" % (program, configSetting))

        else:
            logging.debug("Chewbacca config file not found.  Using defaults.")



    def initalizeCommands(self):
        """
        Provides the static definition of the commandTemplates dictionary.
        """
        # TODO: Mothur doesn't like quoted filenames.  Gross.
        if not ProgramRunner.configsLoaded:
            programPaths = ProgramRunner.programPaths
            ProgramRunner.commandTemplates = {
                "barcode.splitter": "cat \"%s\" | " + programPaths["FASTX"] + 'fastx_barcode_splitter.pl  --bcfile "%s" \
                                    -prefix "%s" --suffix .fastq --bol --mismatches 1',
                "fastx_renamer": programPaths["FASTX"] + "fastx_renamer -n COUNT -i \"%s\" -o \"%s\" -Q 33",
                "pear": programPaths["PEAR"] + " -f \"%s\" -r \"%s\" -o \"%s\" -j %d -m %d",
                "make.contigs": "mothur \'#make.contigs(ffastq=%s, rfastq=%s, bdiffs=1, pdiffs=2, oligos=%s, \
                                    processors=%s)\'",
                "trim.seqs": "mothur \'#trim.seqs(fasta=%s, oligos=%s, maxambig=1, maxhomop=8, \
                                    minlength=300, maxlength=550, bdiffs=1, pdiffs=7)\'",
                "macse_align": "java -jar " + programPaths["MACSE"] + " -prog enrichAlignment  -seq \"%s\" -align \
                                            \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
                                            -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT \"%s\"_NT \
                                            -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",
                "macse_format": "java -jar " + programPaths["MACSE"] + "  -prog exportAlignment -align \"%s\" \
                                            -charForRemainingFS - -gc_def 5 -out_AA \"%s\" -out_NT \"%s\" -statFile \
                                            \"%s\"",
                "trimomatic": "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar SE -%s \"%s\" \"%s\" \
                                            SLIDINGWINDOW:%d:%d MINLEN:%d",
                "chmimera.uchime": "mothur \'#chimera.uchime(fasta=%s, %s)\'",
                "make.fastq": "mothur \'#make.fastq(fasta=%s,qfile=%s)\'",
                "make.fasta": "mothur \'#fastq.info(fastq=%s)\'",
                "remove.seqs": "mothur \'#remove.seqs(accnos=%s, %s)\'",
                "screen.seqs": "mothur \'#screen.seqs(fasta=%s, %s)\'",
                "align.seqs": "mothur \'#align.seqs(candidate=%s, template=%s, flip=t)\'",
                "unique.seqs": "mothur \'#unique.seqs(fasta=%s)\'",
                "cluster-swarm": "",
            }
