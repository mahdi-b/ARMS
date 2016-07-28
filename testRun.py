import subprocess
import glob
import os
rslt = "rslt"
base = "python chewbacca.py "

b = "/home/greg/ARMS/testARMS/testData/barcodes.txt "
f = "/home/greg/ARMS/testARMS/testData/20K_R1.fq "
r = "/home/greg/ARMS/testARMS/testData/20K_R2.fq "
g = "/home/greg/ARMS/testARMS/testData/oligos.txt "
o = "/home/greg/ARMS/src/ARMS/rslt "
odir = o[:-1]




assemble = base + "assemble " + "-n %s -f %s -r %s -o %s -t %d -m %d" % ("test", f,r,o,1,550)
rename = base + "serialize " + "-f %s -r %s -o %s" % (f,r,o)
demux = base + "demux " + "-i %s -b %s -o %s" % ("rslt/test.assembled.fastq",b,o)
makeFasta = [base + "makeFasta" + " -i %s" % x for x in glob.glob(os.path.expanduser("/home/greg/ARMS/src/ARMS/rslt/splitOut*.fastq"))]
trim = [base + "trim " + "-i %s -o %s -g %s" % (x, odir, g) for x in glob.glob(os.path.expanduser("/home/greg/ARMS/src/ARMS/rslt/splitOut*.fasta"))]
tests = [assemble,demux]

try:
    remove = "rm -r rslt"
    #subprocess.check_output(remove, shell=True)
except:
    pass

for test in tests:
    print test
    #subprocess.check_output(test,shell=True)

for test in makeFasta:
    print test
    subprocess.check_output(test, shell=True)

for test in trim:
    print test
    subprocess.check_output(test, shell=True)

