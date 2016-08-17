import ConfigParser
import logging
import os
import subprocess
from Helpers import printVerbose, helpValidate
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
    dry_run = False
    programPaths = {
        "FASTX"  : os.path.expanduser("/usr/bin/fastx_barcode_splitter.pl"),
        "PEAR"   : os.path.expanduser("~/ARMS/programs/pear/pear-0.9.5-bin-64"),
        "USEARCH": os.path.expanduser("~/ARMS/programs/usearch/usearch7.0.1090"),
        "MOTHUR" : os.path.expanduser("~/ARMS/programs/mothur/mothur"),
        "FLEXBAR": os.path.expanduser("~/ARMS/programs/flexbar/flexbar"),
        "VSEARCH": os.path.expanduser("~/ARMS/programs/vsearch/vsearch"),
        "SWARM"  : os.path.expanduser("~/ARMS/programs/swarm/swarm-2.1.9-linux-x86_64")
    }

    def __init__(self, program, params, conditions=None, stdin="", stdout="", stderr=""):
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
        :return             A list of output files to be moved to args.outdir.
        """
        if conditions is None:
            conditions = {}
        if not ProgramRunner.configsLoaded:
            ProgramRunner.loadConfigs(self)
            ProgramRunner.initalizeCommands(self)
            ProgramRunner.configsLoaded = True
        logging.error(params)
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
        #return True
        return helpValidate(conditions)



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


    def run(self):
        """Validates conditions (or prints them for a dry run), sanitizes parameters, and then executes the command.

        :return: None
        """
        try:
            self.validateConditions(self.conditions)
            if printVerbose.VERBOSE or self.dry_run:
                self.dryRun()

            # call and check_call are blocking, Popen is non-blocking
            print( "running " + self.command)
            subprocess.check_call(os.path.expanduser(self.command), shell=True)
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
        config_file_path = ProgramRunner.DEFAULT_CONFIG_FILEPATH
        config_section = "Program Paths"
        if os.path.isfile(config_file_path):
            logging.debug("Loaded Chewbacca config files from %s" % config_file_path)
            config = ConfigParser.RawConfigParser()
            config.read(config_file_path)

            for program in ProgramRunner.programPaths.keys():
                if config.has_option(config_section, program):
                    config_setting = config.get(config_section, program)
                    ProgramRunner.programPaths[program] = os.path.expanduser(config_setting)
                    # logging.debug("Read %s filepath as %s" % (program, config_setting))

        else:
            logging.debug("Chewbacca config file not found.  Using defaults.")



    def initalizeCommands(self):
        """
        Provides the static definition of the commandTemplates dictionary.
        """
        # TODO: Mothur doesn't like quoted filenames.  Gross.
        # TODO: Java jars can't be resolved if string concatenation is used.
        if not ProgramRunner.configsLoaded:
            program_paths = ProgramRunner.programPaths
            ProgramRunner.commandTemplates = {
                "barcode.splitter": "cat \"%s\" | " + program_paths["FASTX"] + '  --bcfile "%s" \
                                            -prefix "%s" --suffix %s --bol --mismatches 1',
                "fastx_renamer": program_paths["FASTX"] + "fastx_renamer -n COUNT -i \"%s\" -o \"%s\" -Q 33",
                "pear": program_paths["PEAR"] + " -f \"%s\" -r \"%s\" -o \"%s\" -v 20 -j %d ",
                "make.contigs": program_paths["MOTHUR"] + " \'#make.contigs(ffastq=%s, rfastq=%s, bdiffs=1, pdiffs=2, \
                                            oligos=%s, processors=%s)\'",
                "trim.seqs": program_paths["MOTHUR"] + " \'#trim.seqs(fasta=%s, oligos=%s, maxambig=1, maxhomop=8, \
                                    minlength=300, maxlength=550, bdiffs=1, pdiffs=7)\'",
                "macse_align": "java -jar ~/ARMS/programs/macse/macse_v1.01b.jar -prog enrichAlignment  -seq \"%s\" \
                                            -align \"%s\" -seq_lr \"%s\" -maxFS_inSeq 0  -maxSTOP_inSeq 0  \
                                            -maxINS_inSeq 0 -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 \
                                            -out_NT \"%s\"_NT -out_AA \"%s\"_AA -seqToAdd_logFile \"%s\"_log.csv",
                "macse_format": "java -jar ~/ARMS/programs/macse/macse_v1.01b.jar  -prog exportAlignment -align \"%s\" \
                                            -charForRemainingFS - -gc_def 5 -out_AA \"%s\" -out_NT \"%s\" -statFile \
                                            \"%s\"",
                "trimmomatic": "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 \
                                            \"%s\" \"%s\" SLIDINGWINDOW:%d:%d MINLEN:%d",
                "chmimera.uchime": program_paths["MOTHUR"] + " \'#chimera.uchime(fasta=%s, %s=%s)\'",
                "make.fastq": program_paths["MOTHUR"] + " \'#make.fastq(fasta=%s,qfile=%s)\'",
                "make.fasta": program_paths["MOTHUR"] + " \'#fastq.info(fastq=%s)\'",
                "remove.seqs": program_paths["MOTHUR"] + " \'#remove.seqs(accnos=%s, %s=%s)\'",
                "screen.seqs": program_paths["MOTHUR"] + " \'#screen.seqs(fasta=%s, %s)\'",
                "flexbar":  program_paths["FLEXBAR"] + " -r \"%s\" -t \"%s\" -ae \"%s\" -a \"%s\"",
                "usearch": program_paths["USEARCH"] + " -derep_fulllength \"%s\" -output \"%s\" -uc \"%s\"",
                "align.seqs": program_paths["MOTHUR"] + " \'#align.seqs(candidate=%s, template=%s, flip=t)\'",
                "vsearch.derep": program_paths["VSEARCH"] + " --derep_fulllength \"%s\" --sizeout --fasta_width 0 "
                                            "--output \"%s\" -uc \"%s\"",
                "swarm": program_paths["SWARM"] + " \"%s\" -o \"%s\" \
                                            -u \"%s\" -w \"%s\"",
                "vsearch.usearch_global": program_paths["VSEARCH"] + " --usearch_global \"%s\" --db \"%s\" --id 0.9 \
                                            --userfields query+target+id+alnlen+qcov --userout \"%s\" --alnout \"%s\"\
                                             %s",

                "unique.seqs": program_paths["MOTHUR"] + " \'#unique.seqs(fasta=%s)\'",
            }
