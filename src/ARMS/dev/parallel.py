from multiprocessing import Pool
from alignReadsProt import *

pool=Pool(2)
senses = [40,20,5]
heuristics = [True, False]
params = [("~/ARMS/data/bold100k.fna bads.fasta", "bads_%d.out" % sens, "~/ARMS/data/bold110k_name_pairs.txt",
          "~/ARMS/data/bold100k.faa", heuristic, "h_%s_rslt.%d.out" % (str(heuristic),sens))
          for sens, heuristic in product(senses, heuristics)]

pool.map_async(run,params).get(9999999)