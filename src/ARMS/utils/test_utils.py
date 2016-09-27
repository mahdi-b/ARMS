from Bio import SeqIO
from Bio.Alphabet import Reduced
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import product

from nose.tools import *

from diffFasta import diff
from src.ARMS.Helpers import makeDirOrdie

test_data = {"A": "a,b,c",
             "B": "b,c,d",
             "C": "",
             "D": ""
             }

test_dir = "testSet"

def parseFastToList(input):
    filetype = input.split('.')[1]
    rslt = []
    for x in SeqIO.parse(open(input, 'r'), filetype):
        rslt.append(x.id)
    return rslt


def makeSeqRecord(id):
    return SeqRecord(Seq("AAA", Reduced.Alphabet), id = id, name = "", description = "",
                     letter_annotations= {"phred_quality":[32,32,32]})

def makeTestFile(file_name, file_type, seq_names):
    with open(file_name, 'w') as out:
        for seq in seq_names:
            if seq != "":
                out.write(makeSeqRecord(seq).format(file_type))
    out.close()
    return file_name


def makeTestSet(file_type):
    makeDirOrdie(test_dir, False)
    test_files = {}
    for entry in test_data.keys():
        sequence_names = test_data[entry].split(",")
        test_files[entry] = makeTestFile("testSet/%s.%s" % (entry, file_type), file_type, sequence_names)
    return test_files

def test_diffFasta():
    def doDiff(X, Y, Z, file_type):
        fasta_file = "%s.%s"
        return diff(X, Y, Z, file_type)

    for file_type in ['fasta','fastq']:
        test_files = makeTestSet(file_type)
        A = test_files["A"]
        B = test_files["B"]
        C = test_files["C"]
        D = test_files["D"]
        rslt = "%s/rslt.%s" % (test_dir,file_type)

        test_set = [A,B,C,D]
        for X,Y in product(test_set, test_set):
            # TODO iterate over test set
            # parse rslt to list
            print (X,Y)
            diff_output = parseFastToList(doDiff(X, Y, rslt, file_type))
            # parse solution to list
            sol = list(set(parseFastToList(X)) - set(parseFastToList(Y)))
            print diff_output
            print sol
            if len(diff_output) == len(sol) and len(sol) == 0 :
                pass
            else:
                assert_equal(diff_output.sort(), sol.sort())


def intersectFasta():
    pass


def test_splitKperFasta():
    pass



def joinFIles():
    pass