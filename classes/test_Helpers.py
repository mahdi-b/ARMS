
import shutil
import StringIO
from Helpers import *
from multiprocessing import Pool
from nose.tools import *

class AttributeObject(object):
    def __init__(self, keys, values):
        for (key, value) in zip(keys, values):
            self.__dict__[key] = value


def mirror(val):
    return val

def test_printVerbose():
    empty_string = ""
    message = "HelloWorld!"
    newline_message = "\n%s" % message
    verbose = [False, True]
    new_line = [False, True]
    rslt = ["", "", message, newline_message]
    i = 0
    for verbosity in verbose:
        printVerbose.VERBOSE = verbosity
        for status in new_line:
            capture = StringIO.StringIO()
            printVerbose(message, newline=status, out=capture)
            assert_equals(capture.getvalue(), rslt[i])
            capture.truncate(0)
            i += 1


def test_helpValidate():
    assert_equals(True, helpValidate({}))
    assert_equals(True, helpValidate({"positive":[1]}))
    assert_raises(Exception, helpValidate({"positive":[0]}))


def test_runProgramRunnerInstance():
    class goodTestRunner:
        def run(self):
            return True
    class badTestRunner:
        def run(self):
            return False
    good = goodTestRunner()
    bad = badTestRunner()

    assert_true(runProgramRunnerInstance(good))
    assert_false(runProgramRunnerInstance(bad))


def test_runPythonInstance():
    vals = [True, False, 3.1415, "test"]
    for val in vals:
        assert_equal(runPythonInstance((mirror,val)), val)


def test_parallel():
    # NOTE: tests for correctness, not parallel execution!!!
    vals = [True, False, 3.1415, "test"]
    rslts = parallel(runPythonInstance, [(mirror, val) for val in vals], Pool(4))
    print rslts
    assert_equal(vals, rslts)


def test_makeDirOrdie():
    test_dir = "1_test_makeDirOrdie_1"
    os.mkdir(test_dir)
    assert_raises(SystemExit, makeDirOrdie,test_dir)
    os.rmdir(test_dir)
    makeDirOrdie(test_dir)
    os.rmdir(test_dir)

def test_makeAuxDir():
    test_dir = "1_test_makeDirOrdie_1"
    aux_dir = "%s_aux" % test_dir
    os.makedirs(aux_dir)
    assert_raises(SystemExit, makeAuxDir,test_dir)
    os.rmdir(aux_dir)
    makeAuxDir(test_dir)
    os.rmdir(aux_dir)

def test_cleanupPool():
    # RUN = 0
    # CLOSE = 1
    # TERMINATE = 2
    pool = Pool(1)
    assert_equal(pool._state, 0)
    assert_raises(SystemExit,cleanupPool,pool)
    assert_equal(pool._state, 2)

def test_strip_ixes():
    dir = "a/b/c"
    name = "tempFile.new.temp.fasta_0"
    ix_string = ""
    ixes=[ "_renamed", "_debarcoded", ".assembled", ".discarded", ".unassembled", "_cleaned", "_derepCount","_derep",
           "_uc", "_splitOut", ".denovo.uchime", "_derepCount", "_uncount", "_counts", "_seeds"]
    for ix in ixes:
        ix_string += ix
    assert_equal(strip_ixes("%s/%s%s%s.ext" % (dir, ix_string, name, ix_string)), name)
    return name


def test_enumerateDir():
    """
    dir = "test_enumerate_Dir"
    files = ["a.txt", "b.fasta", "aa.fasta", "ba.txt", "c"]
    os.mkdir(dir)
    for file_ in files:
        open("%s/%s" % (dir,file_), 'a').close()
    data = [('a.txt',["a.txt"]),
            ('c', ["c"]),
            ('a', []),
            ('a*', ["a.txt", 'aa.fasta']),
            ('*a', ["aa.fasta", "b.fasta"]),
            ('*a*',["a.txt", "aa.fasta", "b.fasta", "ba.txt"] ),
            ('*', files.sort()),
            ('', [])]
    rslts = [enumerateDir(dir,item[0]) for item in data]
    sols = [item[1] for item in data]
    for i in range(len(sols)):
        for j in range(len(sols[i])):
            assert_equal()
        assert_equal()
    shutil.rmtree(dir)
    """
def test_mothur_buildOptionString():
    # Mothur option/attribute names (one filter and one update)
    filterName = "maxlength"
    updateName = "groups"
    # dummy values
    filterProp = 5
    updateProp = "a.group"
    # Test objects with various combinations of filter/update attributes
    emptyObj = AttributeObject([], [])
    filterObj = AttributeObject([filterName], [filterProp])
    updateObj = AttributeObject([updateName], [updateProp])
    completeObj = AttributeObject([filterName, updateName], [filterProp, updateProp])

    for object in [emptyObj, filterObj, updateObj, completeObj]:
        for filter in [True, False]:
            for update in [True, False]:
                try:
                    params = (object, filter, update)
                    attrDict = object.__dict__
                    # print(params)
                    testString = mothur_buildOptionString(*params)
                    print testString
                    # TODO actually check the testString for accuracy
                # if an update file was required, but the object didn't supply one
                except MissingMothurFileException:
                    if not updateName in attrDict:
                        assert True
                    else:
                        print("Unexpected MissingMothurFileException.")
                        assert False

                # if a filter was required but the object didn't supply one
                except MissingMothurFilterException:
                    if not filterName in attrDict:
                        assert True
                    else:
                        print("Unexpected MissingMothurFilterException.")
                        assert False

                # Something went wrong..
                except:
                    print sys.exc_info()
                    assert False

def test_buildCSL():
    testNames = ["a", "b", "c", "d"]
    x = AttributeObject(["a","b","c"],[0,1,None])
    xString = buildCSL(x, testNames)
    assert_equals(xString, "a=0, b=1")
    y = AttributeObject([],[])
    yString = buildCSL(y, testNames)
    assert_equals(yString, "")

test_buildCSL()
test_mothur_buildOptionString()