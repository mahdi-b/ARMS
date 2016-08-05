import sys

biocodeTax = {}


# sys.argv[1] = BiocodePASSED_SAP_tax_info.txt
# sys.argv[2] = vsearch output file with 5 fields


for line in open(sys.argv[1], 'r'):
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


for line in open(sys.argv[2], 'r'):
    data = line.split()
    if biocodeTax.has_key(data[1]):
        data.append(biocodeTax[data[1]])
        print "\t".join(data)
        
    else:
        sys.stderr.write("***********************    Id %s not found ****************************" % data[1])




    
