

class ChewbaccaProgram:
    """Represents a way of completing a ChewbaccaCommand using a particular program (with some combination of external
        programs, and custom Chewbacca scripts.  ChewbaccaPrograms are responsbile for carrying out their associated
        ChewbaccaCommand.  The name attribute denotes the primary external program used, e.g. "vsearch" if VSearch is
        used, or "chewbacca" if a custom script is used.
    """
    # The name of the primary external program used, or 'chewbacca' if a custom script is used.
    name = ""
    # An argparse object
    args = None

    def __init__(self, args_):
        self.args = args_

    def execute_program(self):
        """Completes the representative ChewbaccaCommand using the specified program.
        """
        pass
