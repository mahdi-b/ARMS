import subprocess
from classes.Helpers import makeDirOrdie
from os.path import isdir
from Bio import SeqIO

def cleanup_files(rslt_file):
    cmd = "rm -r %s*" % rslt_file
    try:
        subprocess.check_call(cmd, shell=True)
    except:
        pass


def assert_outdir(outdir):
    return isdir(outdir)


def assert_auxdir(outdir):
    return isdir("%s_aux" % outdir)

def fasta_to_list(path, ext):
    return [x for x in SeqIO.parse(open(path,'r'), ext)]
