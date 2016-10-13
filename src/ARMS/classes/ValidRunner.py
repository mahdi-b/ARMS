from classes import Validator

class ValidRunner(object):
    """A thing that can be validated and run."""

    run_dry = False

    def __init__(self, conditions_=None):
        self.conditions = conditions_
        self.command = ""

    def run(self):
        return


    def dry_run(self, dry_validate=True):
        """Prints the commands that would be run in an actual run, without actually executing them.

        :param dry_validate: If True, print the validation statements that would be executed, but don't actually \
                                execute. If False, perform normal validation routine, exiting on invalid conditions.
        :return: The fully-formatted command string that would be executed.
        """
        return self.command


    def validate_conditions(self):
        """Validates this program's conditions for execution.

        :return:            True if validation is successful, and raising an exception in <Validator.function> if not.
        """
        for condition in self.conditions.items():
            getattr(Validator, condition[0])(condition[1])
        return True


    def dry_validate_conditions(self):
        """Prints validation procedures without actually executing them.
        """
        for condition in self.conditions.items():
            print "\t\tvalidating that %s, %s" % (str(condition[1]), condition[0])


