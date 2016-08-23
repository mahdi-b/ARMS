from chewbaccaFunctions import *
from classes.Helpers import *
from classes.ProgramRunner import ProgramRunner

"""Tests the following chewbaccaFunctions:
        assemble_pear
        splitOnBarcodes
        renameSequences
        trim_flexbar
        trimmomatic
        dereplicate
        align_mothur
        partition
        merge
        ungapFasta
        cluster
        queryBiocode
        queryNCBI
        prescreen
        minhash
"""

def rmPath(path):
    try:
        os.rmdir(path)
    except:
        pass


def assertPathExists(path):
    if not (os.path.isdir(path) or os.path.isfile(path)):
        raise Exception("Path does not exists!")
        exit()


def test_assemble_pear():


def splitOnBarcodes():
    pass


def renameSequences():
    pass


def trim_flexbar():
    pass


def trimmomatic():
    pass


def dereplicate():
    pass


def align_mothur():
    pass


def partition():
    pass


def merge():
    pass


def ungapFasta():
    pass


def cluster():
    pass


def queryBiocode():
    pass


def queryNCBI():
    pass


def prescreen():
    pass


def minhash():
    pass

