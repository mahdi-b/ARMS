import ConfigParser
import logging
import os
import subprocess
from enum import Enum
from Helpers import printVerbose, helpValidate, debugPrint


# Names of available programs
class ProgramRunnerPrograms(Enum):
    FASTX = "FASTX"
    FLEXBAR = "FLEXBAR"
    PEAR = "PEAR"
    SWARM = "SWARM"
    TRIMMOMATIC = "TRIMMOMATIC"
    USEARCH = "USEARCH"
    VSEARCH = "VSEARCH"


# Names of available commands
class ProgramRunnerCommands(Enum):
    ALIGN_VSEARCH = "ALIGN_VSEARCH"
    ASSEMBLE_PEAR = "ASSEMBLE_PEAR"
    CLEAN_TRIMMOMATIC = "CLEAN_TRIMMOMATIC"
    CLUSTER_SWARM = "CLUSTER_SWARM"
    DEMUX_FASTX = "DEMUX_FASTX"
    DEREP_USEARCH = "DEREP_USEARCH"
    DEREP_VSEARCH = "DEREP_VSEARCH"
    TRIM_FLEXBAR = "TRIM_FLEXBAR"
    TEST_ECHO = "TEST_ECHO"


class ProgramRunner(object):
    """A class to interact with external command line programs and internal python functions.  The class contains a \
    dictionary of formatted command strings for external programs.  The class supports validation, sanitization, and \
    debugging of user-supplied parameters.

    Attributes:
        DEFAULT_CONFIG_FILEPATH     The filepath to a config file containing user-specified overrides (such as file
                                        paths to exectuables.
        configsLoaded               Specifies that all config variables have been loaded, and that the file should not
                                        be read again durring the execution of this instance.
        program_paths               Specifies the default location of executibles used in the commands dictionary.
        commandTemplates            A dictionary mapping chewbacca commands to un-paramaterized command line strings.
                                        Used as templates for commands.
    """
    DEFAULT_CONFIG_FILEPATH = "chewbacca.cfg"
    configsLoaded = False
    dry_run = False
    # Actual commands
    commandTemplates = {}

    # Map of Programs to their executables
    program_paths = {
        ProgramRunnerPrograms.FASTX: os.path.expanduser("/usr/bin/fastx_barcode_splitter.pl"),
        ProgramRunnerPrograms.FLEXBAR: os.path.expanduser("~/ARMS/programs/flexbar/flexbar"),
        ProgramRunnerPrograms.PEAR: os.path.expanduser("~/ARMS/programs/pear/pear-0.9.5-bin-64"),
        ProgramRunnerPrograms.SWARM: os.path.expanduser("~/ARMS/programs/swarm/swarm-2.1.9-linux-x86_64"),
        ProgramRunnerPrograms.USEARCH: os.path.expanduser("~/ARMS/programs/usearch/usearch7.0.1090"),
        ProgramRunnerPrograms.VSEARCH: os.path.expanduser("~/ARMS/programs/vsearch/vsearch")
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
        ProgramRunner.loadConfigs(self)
        ProgramRunner.initalizeCommands(self)
        ProgramRunner.configsLoaded = True
        self.program = program
        print self.commandTemplates
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
        for condition in conditions.iteritems():
            print "\t\tvalidating that %s, %s" % (str(condition[1]), condition[0])


    def run(self):
        """Validates conditions (or prints them for a dry run) and then executes the command.

        :return: None
        """

        self.validateConditions(self.conditions)
        if printVerbose.VERBOSE or self.dry_run:
            self.dryRun()

        # call and check_call are blocking, Popen is non-blocking
        debugPrint("Running " + self.command)
        subprocess.check_call(os.path.expanduser(self.command), shell=True)


    def dryRun(self, dryValidate=True):
        """Prints the validation procedures that would be performed, and the commands that would be run in an actual
            run, without actually executing them.

        :return:
        """
        if dryValidate:
            self.dryValidateConditions(self.conditions)
        else:
            self.validateConditions(self.conditions)
        return self.command


    def loadConfigs(self):
        """Loads configuration settings from the chewbacca configuration file (located in
            ProgramRunner.DEFAULT_CONFIG_FILEPATH).  If the file is not found, default values are used. Overwrites the
            default values of any found entries.  Currently loads the [Program  Paths] section, updating the
            program_paths dictionary.
        """
        config_file_path = ProgramRunner.DEFAULT_CONFIG_FILEPATH
        config_section = "Program Paths"
        if os.path.isfile(config_file_path):
            logging.debug("Loaded Chewbacca config files from %s" % config_file_path)
            config = ConfigParser.RawConfigParser()
            config.read(config_file_path)
            for program_name, program_enum in ProgramRunnerPrograms.__members__.iteritems():
                if config.has_option(config_section, program_name):
                    config_setting = config.get(config_section, program_name)
                    ProgramRunner.program_paths[program_enum] = os.path.expanduser(config_setting)
                    # logging.debug("Read %s filepath as %s" % (program, config_setting))

        else:
            logging.debug("Chewbacca config file not found.  Using defaults.")


    def initalizeCommands(self):
        """Provides the static definition of the commandTemplates dictionary.
        """
        # NOTE: Mothur doesn't like quoted filenames.  Gross.
        # NOTE: Java jars can't be resolved if string concatenation is used.

        self.commandTemplates = {
            ProgramRunnerCommands.DEMUX_FASTX: "cat %s | " + self.program_paths[ProgramRunnerPrograms.FASTX] +
                                        '  --bcfile "%s" -prefix "%s" --suffix %s --bol --mismatches 1',
            ProgramRunnerCommands.TEST_ECHO: "echo %s",
            ProgramRunnerCommands.TRIM_FLEXBAR: self.program_paths[ProgramRunnerPrograms.FLEXBAR] +
                                        " -r %s -t %s -ae %s -a %s -u %d",
            ProgramRunnerCommands.ASSEMBLE_PEAR: self.program_paths[ProgramRunnerPrograms.PEAR] +
                                        " -f %s -r %s -o %s -v 20 -j %d ",
            ProgramRunnerCommands.CLUSTER_SWARM: self.program_paths[ProgramRunnerPrograms.SWARM] +
                                        " %s -o %s -u %s -w %s",
            ProgramRunnerCommands.CLEAN_TRIMMOMATIC:
                                        "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar \
                                        SE -phred33 %s %s SLIDINGWINDOW:%d:%d MINLEN:%d",
            ProgramRunnerCommands.DEREP_USEARCH: self.program_paths[ProgramRunnerPrograms.USEARCH] +
                                        " -derep_fulllength %s -output %s -uc %s",
            ProgramRunnerCommands.DEREP_VSEARCH: self.program_paths[ProgramRunnerPrograms.VSEARCH] +
                                        " --derep_fulllength %s --sizeout --fasta_width 0 --output %s -uc %s",
            ProgramRunnerCommands.ALIGN_VSEARCH: self.program_paths[ProgramRunnerPrograms.VSEARCH] +
                                        " --usearch_global %s --db %s --id 0.9 --userfields \
                                        query+target+id+alnlen+qcov --userout %s --alnout %s %s"
        }

