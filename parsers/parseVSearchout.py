import sys

biocodeTax = {}

# sys.argv[1] = BiocodePASSED_SAP_tax_info.txt
# sys.argv[2] = vsearch output file with 5 fields
# sys.argv[3] = min similarity
# sys.argv[4] = min coverage


def parseVSearchout(taxInfo, vsearch_outfile, output_file, min_simmilarity, min_coverage ):
    for line in open(taxInfo, 'r'):
        line = line.rstrip()
        biocodeId = line.split("\t")[0]
        taxonomy = []
        for val in line.split("\t")[1].split(", "):
            data = val.split(": ")
            if len(data) > 1:
                taxonomy.append(data[1])
            else:
                taxonomy.append("")
        biocodeTax[biocodeId] = ":".join(taxonomy)

    min_simm = float(min_simmilarity)
    min_coverage = float(min_coverage)
    with open(output_file, 'w') as out:
        for line in open(vsearch_outfile, 'r'):
            data = line.split()
            # TODO ask mahdi if we need this line
            # the difference in the two files ============================================
            if float(data[2]) < min_simm or float(data[4]) < min_coverage:
                continue
            # end diff ===============================
            if biocodeTax.has_key(data[1]):
                data.append(biocodeTax[data[1]])
                out.write( "\t".join(data))
                out.write("\n")
            else:
                err = "***********************    Id %s not found ****************************" % data[1]
                sys.stderr.write(err)
                out.write(err)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: Biocode_tax_info   vsearch_5_field_output   outfile    min_similarity   min_coverage"
        exit()
    parseVSearchout(sys.argv[1], sys.argv[2])
