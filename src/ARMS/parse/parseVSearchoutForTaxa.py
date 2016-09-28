import sqlite3
import sys

from ete2 import NCBITaxa

from classes.Helpers import printVerbose


# sys.argv[1] = BiocodePASSED_SAP_tax_info.txt
# sys.argv[2] = vsearch output file with 5 fields
# sys.argv[3] = min similarity
# sys.argv[4] = min coverage


def printErrorMissingID(out, ID):
    """Simple error printing for an unresolved/missing reference ID.

    :param out: Outfile handler.
    :param ID: String ID of the unresolved/missing sequence.
    """
    err = "***********************    Id %s not found ****************************\n\
            Missing Taxonomic info for a reference sequence.\n\
            Check that you have provided the correct mapping/reference files.." % ID
    sys.stderr.write(err)
    out.write(err)


def buildTaxaDict(tax_map_file):
    """Takes in a two column tabular file mapping BIOCODE sequence names to semi-colon-delimited taxonomic identifier
        strings and converts it to a dict.

    :param tax_map_file: A two column tabular file mapping BIOCODE sequence names to taxonomic identifier strings.
    :return: A dict where BIOCODE sequence names are keys and map to taxonomic identifier strings.
    """
    biocodeTax = {}
    for line in open(tax_map_file, 'r'):
        line = line.rstrip()
        data = line.split("\t")
        biocodeId = data[0]
        taxonomy = data[1]
        biocodeTax[biocodeId] = taxonomy
    return biocodeTax

# TODO these two be merged
def parseVSearchOutputAgainstFasta(vsearch_outfile, taxInfo, output_file, min_simmilarity, min_coverage):
    """Resolves vsearch matches in a vsearch output file to the taxonomic name taken from BIOCODE.
        Takes in a vsearch output file from usearch__global, parses the result for good matches, and
        writes an output file mapping sequence name to taxa name.

    :param vsearch_out: An output file from vsearch's usearch__global program.
    :param taxInfo: A two column tabular file mapping BIOCODE sequence names to taxonomic identifier strings.
    :param output_file: Where to write the resulting file that maps sequence ID to taxanomic name.
    :param min_coverage: The minimum coverage for an acceptable vsearch match.
    :param min_similarity: The minimum simmilarity for an acceptable vsearch match.
    """
    printVerbose("Parsing Vsearch Output")
    min_simm = float(min_simmilarity)
    min_coverage = float(min_coverage)

    biocodeTax = buildTaxaDict(taxInfo)
    with open(output_file, 'w') as out:
        for line in open(vsearch_outfile, 'r'):
            data = line.split()
            if float(data[2]) < min_simm or float(data[4]) < min_coverage:

                if biocodeTax.has_key(data[1]):
                    print "Found %s as %s" % (data[1], biocodeTax[data[1]])
                    data.append(biocodeTax[data[1]])

                    out.write( "\t".join(data))
                    out.write("\n")
                else:
                    printErrorMissingID(out, data[1])


def parseVSearchOutputAgainstNCBI(vsearch_out, database, output_file, min_coverage, min_similarity):
    """Resolves vsearch matches in a vsearch output file to the taxonomic name taken from BOLD.
        Takes in a vsearch output file from usearch__global, parses the result for good matches, and
        writes an output file mapping sequence name to taxa name.

    :param vsearch_out: An output file from vsearch's usearch__global program.
    :param database: The database used as part of the vsearch usearch__global operation.
    :param output_file: Where to write the resulting file that maps sequence ID to taxanomic name.
    :param min_coverage: The minimum coverage for an acceptable vsearch match.
    :param min_similarity: The minimum simmilarity for an acceptable vsearch match.
    """
    min_simm = float(min_similarity)
    min_coverage = float(min_coverage)
    ncbi = NCBITaxa()
    conn = sqlite3.connect(database)
    c = conn.cursor()

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


    with open(output_file, 'w') as out:
        for line in open(vsearch_out, 'r'):
            data = line.split()

            if float(data[4]) < min_coverage or float(data[2]) < min_simm:
                hit = c.execute(query % data[1]).fetchone()
                if hit:
                    taxonomy = getTaxFromId(hit[0])
                    data.append(taxonomy)
                    print "\t".join(data)
                    out.write("\t".join(data))
                    out.write("\n")
                else:
                    printErrorMissingID(out, data[1])
