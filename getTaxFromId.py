import sys
from ete2 import NCBITaxa
import sqlite3

# sys.argv[1] = $the vsearch output results
# sys.argv[2] = $the ncbi.db database
# sys.argv[3] = min coverage
# sys.argv[4] = min similarity
# the dabase is required to parse the taxonomy from a gi_id

#80972551

ncbi = NCBITaxa()
conn = sqlite3.connect(sys.argv[2])
query = "select taxid from gi_taxid where gi=%s"

def getTaxFromId(taxId, taxonomy=["species", "genus", 'family', 'order', 'class', 'phylum']):
    myTaxonomy = dict([(a,"") for a in taxonomy])
    taxId = int(taxId)
    for lin in ncbi.get_lineage(taxId):

        rank = ncbi.get_rank([lin]).values()[0]
        if rank in taxonomy:
            val = ncbi.get_taxid_translator([lin]).values()[0]
            myTaxonomy[rank] = val

    return ":".join([myTaxonomy[x] for x in taxonomy[::-1]])

c = conn.cursor()

for line in open(sys.argv[1], 'r'):
    data = line.split()

    if float(data[2]) < float(sys.argv[3]) or float(data[4]) < float(sys.argv[4]):
        continue
    hit = c.execute(query % data[1]).fetchone()
    if hit:
        taxonomy = getTaxFromId(hit[0])
        data.append(taxonomy)
        print "\t".join(data)
    else:
        sys.stderr.write("********************    Id %s not found ************************" % data[1])
    






