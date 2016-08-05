import glob
import shutil
import subprocess
import sys
import os
#command  = 
#command = os.path.expanduser(command)

mothur_remove_seqs = ["mothur", "#remove.seqs(accnos=%s, fasta=%s)" % (sys.argv[2], sys.argv[1])]


for file in glob.glob("*.accnos"):
    if len(open(file).readlines()) > 0:
        # run mothur to remove the file
        with open("/tmp/chewbacca.errors", 'a') as fnull:
            print "I am here %s" % mothur_remove_seqs
            subprocess.call(mothur_remove_seqs, stderr=fnull, stdout=fnull)

    else:
        # copy seeds.fasta to seeds.pick.fasta
        shutil.copyfile(sys.argv[1], ".".join(sys.argv[1].split(".")[0:-1])+".pick.fasta")
