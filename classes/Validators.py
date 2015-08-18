import os


def exists(files):
    for file in files:
        if not os.path.exists(file):
            raise Exception("File %s does not exist, executions halted" % file)
    return True
        
def installed(progName, path):
    pass
