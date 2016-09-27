import sys

from classes.Helpers import clip_count


def parseCROPoutToGroups(crop_out_file, output_groups_file):
    """Parses a CROP output file to a groups file.  Crop files are pretty close to the groups format, we just need to
        replace commas with spaces, and clip the counts from child names.
        e.g.
        Crop line:
            BALI4606_0_ID2033_1	BALI4606_0_ID2033_1,BALI4606_0_ID1668_1,BALI4606_0_ID2079_1

        groups line:
            BALI4606_0_ID2033   BALI4606_0_ID2033 BALI4606_0_ID1668 BALI4606_0_ID2079

    :param crop_out_file: The file path to the input crop output file.
    :param output_groups_file: The file path to write the output groups file to.
    :return: The filepath to the output groups file
    """

    with open(output_groups_file, 'w') as out:
        i = 0
        output = ""
        for line in open(crop_out_file, 'r'):
            i += 1
            if i % 100000 == 0:
                out.write(output)
                output = ""
            name, children= line.split("\t")
            seqs = [clip_count(seq) for seq in children.split(",")]
            output += "%s\t%s\n" % (clip_count(name), " ".join(seqs))
        out.write(output)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: crop_out_file  output_groups_file"
    else:
        parseCROPoutToGroups(*sys.argv[1:3])
