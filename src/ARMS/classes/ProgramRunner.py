import ConfigParser
import logging
import os
import subprocess

from enum import Enum

from classes.Helpers import printVerbose, helpValidate, debugPrint


# Names of available programs
class ProgramRunnerPrograms(Enum):
    """An Enum class listing all known external programs used by chewbacca
    """
    CROP = "CROP"
    FASTX = "FASTX"
    FLEXBAR = "FLEXBAR"
    PEAR = "PEAR"
    SPADES = "SPADES"
    SWARM = "SWARM"
    TRIMMOMATIC = "TRIMMOMATIC"
    VSEARCH = "VSEARCH"
    MACSE = "MACSE"
    JAVA = "JAVA"


# Names of available commands
class ProgramRunnerCommands(Enum):
    """An Enum class listing all known external commandline operations."""
    ALIGN_VSEARCH = "ALIGN_VSEARCH"
    ASSEMBLE_PEAR = "ASSEMBLE_PEAR"
    CLEAN_TRIMMOMATIC = "CLEAN_TRIMMOMATIC"
    CLUSTER_CROP = "CLUSTER_CROP"
    CLUSTER_SWARM = "CLUSTER_SWARM"
    CLUSTER_VSEARCH = "CLUSTER_VSEARCH"
    DEMUX_FASTX = "DEMUX_FASTX"
    DEREP_VSEARCH = "DEREP_VSEARCH"
    PRECLEAN_SPADES = "PRECLEAN_SPADES"
    TRIM_FLEXBAR = "TRIM_FLEXBAR"
    TEST_ECHO = "TEST_ECHO"
    MACSE_ALIGN = "MACSE_ALIGN"
    MACSE_FORMAT = "MACSE_FORMAT"


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
    run_dry = False
    # Actual commands
    commandTemplates = {}

    # Map of Programs to their executables
    program_paths = {
        ProgramRunnerPrograms.CROP: os.path.expanduser("~/ARMS/programs/crop/crop"),
        ProgramRunnerPrograms.FASTX: os.path.expanduser("/usr/bin/fastx_barcode_splitter.pl"),
        ProgramRunnerPrograms.FLEXBAR: os.path.expanduser("~/ARMS/programs/flexbar/flexbar"),
        ProgramRunnerPrograms.PEAR: os.path.expanduser("~/ARMS/programs/pear/pear-0.9.5-bin-64"),
        ProgramRunnerPrograms.SPADES: os.path.expanduser("~/ARMS/programs/spades/bin/spades.py"),
        ProgramRunnerPrograms.SWARM: os.path.expanduser("~/ARMS/programs/swarm/swarm-2.1.9-linux-x86_64"),
        ProgramRunnerPrograms.VSEARCH: os.path.expanduser("~/ARMS/programs/vsearch/vsearch"),
        ProgramRunnerPrograms.TRIMMOMATIC: os.path.expanduser("~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar"),
        ProgramRunnerPrograms.MACSE: os.path.expanduser("~/ARMS/programs/macse/macse_v1.01b.jar"),
        ProgramRunnerPrograms.JAVA: os.path.expanduser("java -jar")
    }

    def __init__(self, program, params, conditions=None, custom_arg_string="", stdin="", stdout="", stderr=""):
        """Initalizes a ProgramRunner object.  Reads chewbacca.cfg and loads configuration settings, builds the command
        string, and configures stdIO pipes.

        :param program:     A chewbacca command (as a string).  Used as a key to fetch the command string from
                                ProgramRunner.commandTemplates.
        :param params:      A list of parameters for the chosen chewbacca command.
        :param conditions:  A dictionary mapping <Validator>.function names to a list of parameters.  The dictionary
                                key/function name is called on each parameter in the corresponding list of parameters.
                                All parameters must satisfy their <Validator>.function in order for the specified
                                program to execute.
        :param custom_arg_string: A string containing advanced command line parameters to pass to the called program.
                                This string will be appended to the end of the normal commandline arguments.
        :param stdin:       A handle to the stdin pipe for the command.  Defaults to stdin.
        :param stdout:      A handle to the stdout pipe for the command.  Defaults to stdout.
        :param stderr:      A hande to the stderr pipe for the command.  Defaults to stderr.
        :return             A list of output files to be moved to outdir.
        """
        if conditions is None:
            conditions = {}
        ProgramRunner.loadConfigs(self)
        ProgramRunner.initalizeCommands(self)
        ProgramRunner.configsLoaded = True
        self.program = program
        self.command = "%s %s" % (self.commandTemplates[program] % tuple(params), custom_arg_string)
        self.conditions = conditions
        self.extra_args = custom_arg_string
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr

    def validate_conditions(self, conditions):
        """Validates this program's conditions for execution.

        :param conditions: A dictionary mapping <Validator>.function names to a list of parameters.  The dictionary
                            key/function name is called on each parameter in the corresponding list of parameters.
                            All parameters must satisfy their <Validator>.function in order for the specified
                            program to execute.

        :return:            True if validation is successful, and raising an exception in <Validator.function> if not.
        """
        return helpValidate(conditions)

    def dry_validate_conditions(self, conditions):
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

        self.validate_conditions(self.conditions)
        if printVerbose.VERBOSE or self.run_dry:
            self.dry_run()

        # call and check_call are blocking, Popen is non-blocking
        debugPrint("Running " + self.command)
        subprocess.check_call(os.path.expanduser(self.command), shell=True)

    def dry_run(self, dryValidate=True):
        """Prints the validation procedures that would be performed, and the commands that would be run in an actual
            run, without actually executing them.

        :return:
        """
        if dryValidate:
            self.dry_validate_conditions(self.conditions)
        else:
            self.validate_conditions(self.conditions)
        return self.command

    def loadConfigs(self):
        """Loads configuration settings from the chewbacca configuration file (located in
            ProgramRunner.DEFAULT_CONFIG_FILEPATH).  If the file is not found, default values are used. Overwrites the
            default values of any found entries.  Currently loads the [Program  Paths] section, updating the
            program_paths dictionary.
        """
        bad_configs = False
        config_file_path = ProgramRunner.DEFAULT_CONFIG_FILEPATH
        config_section = "Program Paths"
        if os.path.isfile(config_file_path):
            logging.debug("Loaded Chewbacca config files from %s" % config_file_path)
            config = ConfigParser.RawConfigParser()
            config.read(config_file_path)
            for program_name, program_enum in ProgramRunnerPrograms.__members__.iteritems():
                if config.has_option(config_section, program_name):
                    config_setting = config.get(config_section, program_name)
                    if os.path.isfile(os.path.expanduser(config_setting)):
                        ProgramRunner.program_paths[program_enum] = os.path.expanduser(config_setting)
                    else:
                        print "ERROR! The provided config path does not exist: %s = %s" % (program_name, config_setting)
                        bad_configs = True
                    # logging.debug("Read %s filepath as %s" % (program, config_setting))
            if bad_configs:
                print "Errors found in configs file.  Please check config settings."
                exit()
        else:
            logging.debug("Chewbacca config file not found.  Using defaults.")

    def initalizeCommands(self):
        """Initalizes the commandTemplates dictionary.
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
            ProgramRunnerCommands.CLUSTER_CROP: self.program_paths[ProgramRunnerPrograms.CROP] +
                                                 " -i %s -o %s -z %d -%s -e %d -m %d -r %d %s",
            ProgramRunnerCommands.CLUSTER_SWARM: self.program_paths[ProgramRunnerPrograms.SWARM] +
                                        " %s -o %s -u %s -w %s",
            ProgramRunnerCommands.CLUSTER_VSEARCH: self.program_paths[ProgramRunnerPrograms.VSEARCH] +
                                        " --cluster_size %s -id %f --centroids %s --uc %s",
            ProgramRunnerCommands.CLEAN_TRIMMOMATIC: self.program_paths[ProgramRunnerPrograms.JAVA] +
                                        " " + self.program_paths[ProgramRunnerPrograms.TRIMMOMATIC] +
                                        " " + "SE -phred33 %s %s SLIDINGWINDOW:%d:%d MINLEN:%d",
            ProgramRunnerCommands.DEREP_VSEARCH: self.program_paths[ProgramRunnerPrograms.VSEARCH] +
                                        " --threads %d --derep_fulllength %s --sizeout --fasta_width 0 --output %s \
                                        -uc %s",
            ProgramRunnerCommands.ALIGN_VSEARCH: self.program_paths[ProgramRunnerPrograms.VSEARCH] +
                                        " --threads %d --usearch_global %s --db %s --id 0.9 --userfields \
                                        query+target+id+alnlen+qcov --userout %s --alnout %s %s",
            ProgramRunnerCommands.PRECLEAN_SPADES: "python " + self.program_paths[ProgramRunnerPrograms.SPADES] +
                                         " --only-error-correction --disable-gzip-output -1 %s -2 %s -o %s -t %d",
            ProgramRunnerCommands.MACSE_ALIGN:  self.program_paths[ProgramRunnerPrograms.JAVA] +
                                        " " + self.program_paths[ProgramRunnerPrograms.MACSE] +
                                        " -prog enrichAlignment  -seq %s -align \
                                        %s -seq_lr %s -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
                                        -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT %s_NT \
                                        -out_AA %s_AA -seqToAdd_logFile %s_log.csv",
            ProgramRunnerCommands.MACSE_FORMAT: self.program_paths[ProgramRunnerPrograms.JAVA] +
                                        " " + self.program_paths[ProgramRunnerPrograms.MACSE] +
                                        " " + " -prog exportAlignment -align %s \
                                        -charForRemainingFS - -gc_def 5 -out_AA %s -out_NT %s -statFile %s",
}
