import subprocess
import glob
import os
rslt = "rslt"
base = "python chewbacca.py "
currstep = 0
b = "/home/greg/ARMS/testARMS/testData/barcodes.txt "
f = "/home/greg/ARMS/testARMS/testData/20K_R1.fq "
r = "/home/greg/ARMS/testARMS/testData/20K_R2.fq "
g = "/home/greg/ARMS/testARMS/testData/oligos.txt "

assembler_prog = "pear"

def lastStep():
    global currstep
    return currstep - 1

def nextStep():
    global currstep
    currstep = currstep + 1
    return "step%d" % currstep

def reset(rslt):
    try:
        remove = "rm -r rslt"
        subprocess.check_output(remove, shell=True)
    except:
        pass


assemble = base + "assemble " + "-p %s -n %s -f %s -r %s -o %s -t %d -m %d" % (assembler_prog, "test", f, r, "assemble", 1, 550)
rename = base + "serialize " + "-f %s -r %s -o %s" % (f, r, "rename")

demux = base + "demux " + "-i %s -b %s -o %s" % ("assemble/test.assembled.fastq", b, "demux")

makeFasta = [base + "makeFasta" + " -i %s" % x for x in glob.glob(os.path.expanduser("/home/greg/ARMS/src/ARMS/rslt/splitOut*.fastq"))]

trim = [base + "mothurTrim " + "-i %s -o %s -g %s" % (x, "trimmed", g) for x in glob.glob("2_split/splitOut*.*")]


#reset()
#subprocess.check_output(assemble,shell=True)
#subprocess.check_output(demux, shell=True)
#for t in trim:
#    subprocess.check_output(t, shell=True)

print reset
print assemble
print demux
for t in trim:
	print t


