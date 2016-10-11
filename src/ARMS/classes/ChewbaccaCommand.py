

class ChewbaccaCommand:
    """Represents a portable, executable Command.  ChewbaccaCommands have a list of supported_programs, a
        default_program, and an args object (object with parameters as attributes).  Chewbacca Commands read the
        program attribute from the args object, and call the appropriate program (if supported), falling back to the
        default if no match is found.
    """
    supported_programs = []
    default_program = None
    command_name = ""
    args = None

    def __init__(self, args_):
        self.args = args_

    def get_program(self, program):
        """Given a program name, searches for that Program in the list of supported_programs, and returns an instance
        of that Program initalized with the parameter object.  If not found, returns an instance of the default_program
        initalized with the parameter object.

        :param program: String name of the desired program.  This should match the .name attribute of a ChewbaccaProgram
        :return: An instance of a ChewbaccaProgram, initalized with the argparse object.
        """
        for prog in self.supported_programs:
            if prog.name == program:
                print "Executing command %s with program %s." % (self.command_name, prog.name)
                return prog(self.args)
        print "Command %s has no program '%s'.  Using default program %s." % (self.command_name, program,
                                                                              self.default_program.name)
        return self.default_program(self.args)

    def execute_command(self):
        """Executes the command.

        :return: Results of execution.
        """
        program_name = getattr(self.args, 'program', "")

        return self.get_program(program_name).execute_program()
