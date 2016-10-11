from classes.ValidRunner import ValidRunner
from classes.Helpers import printVerbose

class PythonRunner(ValidRunner):
    def __init__(self, function_, params_, conditions_):
        self.function = function_
        self.params = params_
        self.conditions = conditions_
        self.command = "python <%s>: %s" % (function_.__name__, self.params)

    def run(self):
        printVerbose("Running: %s with %s" % (self.function.__name__, str(self.params)))
        return self.function(*self.params)
