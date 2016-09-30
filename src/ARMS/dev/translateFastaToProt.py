import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from itertools import product

table_list = [5, 9, 4, 2, 1, 13, 14, 6, 3, 10, 11, 12]
orf = range(3)

def parse_gc_table(input_file="uniq_short.txt", delim="\t", gc_col=0, tax_col=1, tax_start=2, tax_end=3,tax_delim=";"):

    gc_lookup = {}
    for line in open(input_file, 'r'):
        data = line.split(delim)
        tax = data[tax_col]
        gc = int(data[gc_col])
        name = ";".join(tax.strip().split(tax_delim)[tax_start:tax_end])
        gc_lookup[name] = gc
    return gc_lookup


def getTaxName(record, depth=3):
    id = str(record.description)
    # lines look like:
    # seq_id<tab>order:data,family:data,genus:data...
    seqID, tax = id.split("\t")
    tax_data = tax.split(",")[0:depth]
    tax_parts = [data.split(":")[1] for data in tax_data]
    tax_name = ";".join(tax_parts)
    return tax_name


def translate_fasta(input_file, gc_file, filetype="fasta"):
    possibilities = defaultdict(int)
    perfect_translations = []
    unknown_translations = []
    table = {}
    prefix = input_file.split(".")[0]
    for i in range(3):
        table.update(parse_gc_table(input_file = gc_file, tax_start=0, tax_end=i+1))
    with open("%s_translated.%s" % (prefix, filetype), 'w') as good_output:
        with open("%s_unknown.%s" % (prefix, filetype), 'w') as bad_output:
            i = 0
            for record in SeqIO.parse(open(input_file, "r"), filetype):
                i += 1
                if i % 10000 == 0:
                    SeqIO.write(perfect_translations, good_output, filetype)
                    perfect_translations = []
                    SeqIO.write(unknown_translations, bad_output,  filetype)
                    unknown_translations = []

                my_tax = getTaxName(record)
                my_gc = 5
                if table.has_key(my_tax):
                    my_gc = table[my_tax]
                possible_translations = translate_with_known_gc(record, [my_gc])

                if len(possible_translations) == 0:
                    possible_translations = translate_with_unknow_gc(record)

                if len(possible_translations) == 1:
                    record.seq = Seq(possible_translations[0])
                    perfect_translations.append(record)
                else:
                    unknown_translations.append(record)
                possibilities[len(possible_translations)] += 1

            SeqIO.write(perfect_translations, good_output, filetype)
            SeqIO.write(unknown_translations, bad_output, filetype)

    print possibilities


def translate_with_known_gc(record, gc_list):
    return try_translate(record, gc_list)


def translate_with_unknow_gc(record):
    return try_translate(record)


def try_translate(record, gc_list = table_list):
    possible_translations = []
    for orf, table in product(range(3), gc_list):
        translation = str(record[orf:].seq.translate(table=table))
        rc_translation = str(record[orf:].reverse_complement().seq.translate(table=table))
        # If no stop codons, we're done
        if translation.count("*") == 0:
            possible_translations .append(translation)
        if rc_translation.count("*") == 0:
            possible_translations.append(rc_translation)
    return possible_translations


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: input_file  gc_table,"
    else:
        translate_fasta(*sys.argv[1:3])
    #bestHits = findBestHits("%s%s" % (data_dir, "/data/test2"))
