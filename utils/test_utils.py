from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import Reduced
from nose.tools import *

test_data = {"A": "a,b,c",
             "B": "b,c,d",
             "C": "",
             "D": ""
             }

test_dir = "testSet"

def parseFastToList(input, filetype):
    return list(SeqIO.parse(open(input,'r'), filetype))

def makeSeqRecord(id):
    return SeqRecord(Seq("AAA", Reduced.Alphabet), id = id, name = "", description = "")

def makeTestFile(file_name, file_type, seq_names):
    with open("%s.%s" % (file_name, file_type), 'w') as out:
        for seq in seq_names:
            out.write(makeSeqRecord(seq).format(file_type))
    out.close()
    return file_name
def makeTestSet(file_type):
    test_files = {}
    for entry in test_data.keys():
        test_files[entry] = makeTestFile("testSet/%s" % entry, file_type, test_data[entry].split(","))
    return test_files

def test_diffFasta():
    for file_type in ['fasta','fastq']:
        test_files = makeTestSet('fasta')
        A = test_files["A"]
        B = test_files["B"]
        C = test_files["C"]
        D = test_files["D"]
        rslt = "%s/rslt.%s" % (test_dir, file_type)

        # TODO iterate over test set
        assert_equal(list(set(parseFastToList(A,B, rslt))) , parseFastToList(rslt, file_type))

    pass


def intersectFasta():
    pass


def test_splitKperFasta():
    pass



def joinFIles():
    pass