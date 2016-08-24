import sqlite3
import sys
from ete2 import NCBITaxa
# sys.argv[1] = $the vsearch output results
# sys.argv[2] = $the ncbi.db database
# sys.argv[3] = min coverage
# sys.argv[4] = min similarity
# the dabase is required to parse the taxonomy from a gi_id

#80972551
def parseVSearchToTaxa(vsearch_out, ncbi_db, output_file, min_coverage, min_similarity):
    ncbi = NCBITaxa()
    conn = sqlite3.connect(ncbi_db)
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
    min_simm = float(min_similarity)
    min_coverage = float(min_coverage)
    with open(output_file, 'w') as out:
        for line in open(vsearch_out, 'r'):
            data = line.split()

            if float(data[4]) < min_coverage or float(data[2]) < min_simm:
                continue
            hit = c.execute(query % data[1]).fetchone()
            if hit:
                taxonomy = getTaxFromId(hit[0])
                data.append(taxonomy)
                print "\t".join(data)
                out.write("\t".join(data))
                out.write("\n")
            else:
                err = "********************    Id %s not found ************************" % data[1]
                sys.stderr.write(err)
                out.write(err + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print "Usage: vsearch_5_field_output   ncbi_db   outfile     min_similarity   min_coverage"
    else:
        parseVSearchToTaxa(*sys.argv[1:6])




