import ConfigParser
import logging
import os
import subprocess
import sys
from enum import Enum
from classes.Helpers import printVerbose, debugPrint
from classes.ValidRunner import ValidRunner

class ProgramRunnerPrograms(Enum):
    """An Enum class listing all known external programs used by chewbacca."""
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


class ProgramRunner(ValidRunner):
    """A class to interact with external command line programs.  The class contains a \
    dictionary of formatted command strings for external programs.  The class supports validation, sanitization, and \
    debugging of user-supplied parameters.

    Attributes:
        self.DEFAULT_CONFIG_FILEPATH:     The filepath to a config file containing user-specified overrides (such as \
                                            file paths to exectuables.
        self.program_paths:               Specifies the default location of executables used in the commandTemplates \
                                            dictionary.
        self.commandTemplates:            A dictionary mapping chewbacca commands to un-paramaterized command line \
                                            strings.  Used as templates for commands.
    """
    DEFAULT_CONFIG_FILEPATH = "chewbacca.cfg"

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

    def __init__(self, program_, params, conditions_={}, custom_arg_string="", dryrun=False):
        """Initalizes a ProgramRunner object.  Reads chewbacca.cfg and loads configuration settings, builds the command
        string, and configures stdIO pipes.

        :param program_:     A ProgramRunnerCommand.  Used as a key to fetch the command string from
                                ProgramRunner.commandTemplates.
        :param params:      A list of parameters for the chosen chewbacca command.
        :param conditions_:  A dictionary mapping <Validator>.function names to a list of parameters.  The dictionary
                                key/function name is called on each parameter in the corresponding list of parameters.
                                All parameters must satisfy their <Validator>.function in order for the specified
                                program to execute.
        :param custom_arg_string: A string containing advanced command line parameters to pass to the called program.
                                This string will be appended to the end of the normal commandline arguments.
        :return             A list of output files to be moved to outdir.
        """
        if conditions_ is None:
            conditions = {}
        ProgramRunner.__load_configs__(self)
        ProgramRunner.__initalize_commands__(self)
        self.program = program_
        self.command = "%s %s" % (self.commandTemplates[program_] % tuple(params), custom_arg_string)
        self.conditions = conditions_
        self.extra_args = custom_arg_string
        self.run_dry = dryrun

    def run(self):
        """Validates conditions (or prints them for a verbose/dry run) and then executes the command.
        """

        # If verbose or dry run, print procedures, otherwise, silence stdout.
        if self.run_dry:
            return self.dry_run()
        else:
            if printVerbose.VERBOSE:
                output = sys.stdout
                self.dry_validate_conditions()
            else:
                output = open(os.devnull, 'w')

            self.validate_conditions()
            debugPrint("Running " + self.command)
            # call and check_call are blocking, Popen is non-blocking
            try:
                rslt = subprocess.check_call(os.path.expanduser(self.command), shell=True, stdout=output)
            except subprocess.CalledProcessError as cpe:
                pass
            return

    def __load_configs__(self):
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
            for program_name, program_enum in ProgramRunnerPrograms.__members__.items():
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

    def __initalize_commands__(self):
        """Initalizes the commandTemplates dictionary.
        """

        self.commandTemplates = {
            ProgramRunnerCommands.DEMUX_FASTX: "cat %s | " + self.program_paths[ProgramRunnerPrograms.FASTX] +
                                                    "  --bcfile %s -prefix %s --suffix %s --bol --mismatches 1",
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
                                                    " SE -phred33 %s %s SLIDINGWINDOW:%d:%d MINLEN:%d",
            ProgramRunnerCommands.DEREP_VSEARCH: self.program_paths[ProgramRunnerPrograms.VSEARCH] +
                                                    " --threads %d --quiet --derep_fulllength %s --sizeout \
                                                    --fasta_width 0 --output %s -uc %s",
            ProgramRunnerCommands.ALIGN_VSEARCH: self.program_paths[ProgramRunnerPrograms.VSEARCH] +
                                                    " --threads %d --usearch_global %s --db %s --id %f --userfields \
                                                    query+target+id+alnlen+qcov+qstrand --userout %s --alnout %s %s",
            ProgramRunnerCommands.PRECLEAN_SPADES: "python " + self.program_paths[ProgramRunnerPrograms.SPADES] +
                                                    " --only-error-correction --disable-gzip-output -1 %s -2 %s -o %s \
                                                    -t %d",
            ProgramRunnerCommands.MACSE_ALIGN: self.program_paths[ProgramRunnerPrograms.JAVA] +
                                                    " " + self.program_paths[ProgramRunnerPrograms.MACSE] +
                                                    " -prog enrichAlignment  -seq %s -align \
                                                    %s -seq_lr %s -maxFS_inSeq 0  -maxSTOP_inSeq 0  -maxINS_inSeq 0 \
                                                    -maxDEL_inSeq 3 -gc_def 5 -fs_lr -10 -stop_lr -10 -out_NT %s_NT \
                                                    -out_AA %s_AA -seqToAdd_logFile %s_log.csv",
            ProgramRunnerCommands.MACSE_FORMAT: self.program_paths[ProgramRunnerPrograms.JAVA] +
                                                    " " + self.program_paths[ProgramRunnerPrograms.MACSE] +
                                                    " -prog exportAlignment -align %s \
                                                    -charForRemainingFS - -gc_def 5 -out_AA %s -out_NT %s -statFile %s",
        }
