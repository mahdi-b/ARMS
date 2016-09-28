

class ChewbaccaCommand:
    supported_programs = []
    default_program = None
    command_name = ""
    args = None

    def __init__(self, args_):
        self.args = args_


    def get_program(self, program):
        for prog in self.supported_programs:
            if prog.name == program:
                print "Executing command %s with program %s." % (self.command_name, prog.name)
                return prog(self.args)
        return self.default_program(self.args)


    def execute_command(self):
        pass