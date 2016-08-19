import sys
import glob
import os
from collections import defaultdict
from classes.Helpers import *
seqTosample=dict()
samples =set()

#python ~/bin/builMatrix.py ../5_dereplicate/ ../9_cluster/uc_parsed.out ../9_cluster/clustering.out matrix.out
def buildMatrix(dereplicated_names_file, pre_cluster_derep_names, post_clustering_names_file, out_file):
    for f in getInputs(dereplicated_names_file):
        sampleName = "_".join(os.path.basename(f).split("_")[1:-4])
        samples.add(sampleName)
        for line in open(f).readlines():
            line = line.rstrip()
            data = line.split("\t")
            seqTosample[data[0]] = sampleName
            if len(data) == 1:
                continue
            for s in data[1].split():
                seqTosample[s] = sampleName

    derepInfoLevel_II = dict()

    for line in open(pre_cluster_derep_names):
        line = line.rstrip()
        data = line.split("\t")
    
        if len(data) > 1:
            derepInfoLevel_II[data[0].split("_")[0]] = [data[0]]+data[1].split()
        else:
            derepInfoLevel_II[data[0].split("_")[0]] = [data[0]]
    
    outFile = open(out_file, 'w')
    outFile.write("\t"+"\t".join(samples)+'\n')
    
    for line in open(post_clustering_names_file):
        data = line.split()
        counts = defaultdict(int)
        seq_to_sample_keys = seqTosample.keys()
        for d in data:
            for e in derepInfoLevel_II[d.split("_")[0]]:
                id = e.split("_")[0]
                if id in seq_to_sample_keys:
                    sample = seqTosample[id]
                    if sample in counts.keys():
                        counts[sample] += int(e.split("_")[1])
    
        outFile.write("%s\t" % data[0].split("_")[0])
        for sample in samples:
            outFile.write("%s\t" % counts[sample])
        outFile.write("\n")
    
    outFile.close()

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Usage: input_dir arg2 arg3 out_file"
        exit()
    buildMatrix(*sys.argv[1:5])