from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import Reduced
from nose.tools import *
from ..classes.Helpers import makeDirOrdie
from diffFasta import diff
test_data = {"A": "a,b,c",
             "B": "b,c,d",
             "C": "",
             "D": ""
             }

test_dir = "testSet"

def parseFastToList(input, filetype):
    rslt = []
    for x in SeqIO.parse(open(input, 'r'), filetype):
        rslt.append(x)
    return rslt


def makeSeqRecord(id):
    return SeqRecord(Seq("AAA", Reduced.Alphabet), id = id, name = "", description = "")

def makeTestFile(file_name, file_type, seq_names):
    with open("%s.%s" % (file_name, file_type), 'w') as out:
        for seq in seq_names:
            out.write(makeSeqRecord(seq).format(file_type))
    out.close()
    return file_name


def makeTestSet(file_type):
    makeDirOrdie(test_dir, False)
    test_files = {}
    for entry in test_data.keys():
        test_files[entry] = makeTestFile("testSet/%s.%s" % entry, file_type, test_data[entry].split(","))
    return test_files

def test_diffFasta():
    def doDiff(X, Y, Z, file_type):
        fasta_file = "%s.%s"
        return diff(fasta_file % (X , file_type), fasta_file % (Y , file_type), fasta_file % (Z , file_type))

    for file_type in ['fasta','fastq']:
        test_files = makeTestSet('fasta')
        A = test_files["A"]
        B = test_files["B"]
        C = test_files["C"]
        D = test_files["D"]
        rslt = "%s/rslt.%s" % (test_dir, file_type)

        # TODO iterate over test set
        # parse rslt to list
        print parseFastToList(doDiff(A,B, rslt, file_type), file_type)
        # parse solution to list
        print set(parseFastToList(A, file_type)) - set(parseFastToList(B, file_type))


    pass


def intersectFasta():
    pass


def test_splitKperFasta():
    pass



def joinFIles():
    pass