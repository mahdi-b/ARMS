

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
import os
from chewbaccaFunctions import *
from classes import Helpers
from nose.tools import *


test_data_dir = os.path.expanduser("~/ARMS/data")
test_outdir = os.path.expanduser("~/ARMS/testARMS/testData/")
emptyFile = getInputs(test_data_dir, "empty.file")


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
    out_dir = "%s/%s" % (test_outdir, "1_assemble")
    test_fn = assemble_pear
    f_read = test_data_dir
    r_read = test_data_dir
    bad_dir = os.path.expanduser("~")
    def call_fn(f, r, outdir):
        rmPath(out_dir)
        args = lambda: 0
        args.input_f = f
        args.input_r = r
        args.name = "BALI"
        args.outdir = outdir
        args.threads = 1
        test_fn(args)

    def test_inputs(f, r, outdir):
        # If missing an input or output, exit
        assert_raises(SystemExit, call_fn, bad_dir, r_read, out_dir)
        assert_raises(SystemExit, call_fn, f_read, bad_dir, out_dir)
        assert_raises(SystemExit, call_fn, f_read, r_read, bad_dir)
        # Valid call
        call_fn(f_read, r_read, out_dir)

    f_read = getInputs(test_data_dir, "*_R1*")[0]
    r_read = getInputs(test_data_dir, "*_R1*")[0]


def splitOnBarcodes():
    test_fn = splitOnBarcodes
    input = test_data_dir
    outdir = test_outdir
    barcodes = emptyFile
    def call_fn(input, barcodes, outdir):
        rmPath(test_outdir)
        args = lambda: 0
        args.input = input
        args.output = outdir
        args.barcodes = barcodes
        test_fn(args)

    # If missing an input or output, exit
    assert_raises(SystemExit, call_fn, "", outdir, barcodes)
    assert_raises(SystemExit, call_fn, input, "", test_outdir)
    assert_raises(SystemExit, call_fn, input, outdir, "")


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

