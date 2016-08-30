import sys
import re


# TODO REMOVE THIS
def flag_frame_shift_seqs(inputFile):
    """Flags sequences with a frameshift in a caln file.

    :param inputFile:   Input file in caln format.  Each line contains some permutation of digits, with I,D,M.
    :return:    Returns a tuple, where the first element is an array of sequences which likely contain a frameshift,
                    and the second element is the total number of lines (and therefore sequences) in teh inputFile.
    """
    current_line = 0
    prefix_indels = 0
    suffix_indels = 0
    with open(inputFile, 'r') as input:
        flagged_seqs = []
        for line in input:
            flagged = False
            rslt = re.search(r'^(\d+)[ID].*M(\d+)[ID]$', line)
            if rslt != None:
                prefix_indels, suffix_indels = rslt.groups(0)
            # sequence_length = re.search(r'^\d*[ID](\d*)M\d*[ID]$', line).groups(0)
            if re.search(r'.*M.*[ID]+.*M', line):
                flagged = True
            """
            elif int(prefix_indels) % 3 != 0:
                flagged = True
            elif int(suffix_indels)  %3 != 0:
                flagged = True
            """
            if flagged:
                flagged_seqs.append(current_line)
            current_line = current_line + 1
    input.close()
    return flagged_seqs, current_line


def scrape_ulog(inputFile, outputFile, target_lines):
    """Scrapes data from an output Vsearch .aln file, into an outputFile.  Entries to be scraped are listed in a list
        of entry numbers in target_lines.

    :param inputFile:       The output .aln file from vsearch.
    :param outputFile:      The output file path to write to.
    :param target_lines:    A list of entry numbers to scrape from the inputFile.
    :return:
    """
    if len(target_lines) == 0:
        print "Nothing to scrape"
        return

    query_name = ""
    target_name = ""
    query_seq = ""
    target_seq = ""
    match_seq = ""
    tmp = ""
    query_start = 0
    target_start = 0
    query_end = 0
    target_end = 0
    match_start = 0
    new_target = True
    new_query = True
    read_match = False
    current_entry = 0
    current_target_index = 0
    entries_to_pick = sorted(target_lines)
    with open(inputFile, 'r') as input:
        with open(outputFile, 'w') as output:
            #if current_entry == entries_to_pick[current_target_index]:
            for line in input:
                """
                print "=============="
                print line
                print "=============="
                """
                # if our current entry count is larger than the last target index, we're done
                if current_entry > entries_to_pick[-1]:
                    print ("Terminated Early @ %dth entry of input file" % current_entry)
                    output.close()
                    input.close()
                    return
                if line == "" or line == "\n":
                    pass
                elif line[:6] == " Query" and current_entry == target_lines[current_target_index]:
                    #print(1)
                    query_name = line[line.find('>') + 1:].rstrip('\n')
                    query_seq = ""
                    match_seq = ""
                    new_query = True

                elif line[:6] == "Target" and current_entry == target_lines[current_target_index]:
                    #print(2)
                    target_name = line[line.find('>') + 1:].rstrip('\n')
                    #print target_name
                    target_seq = ""
                    new_target = True

                elif line[:3] == "Qry" and current_entry == target_lines[current_target_index]:

                    #print(3)
                    query_seq += line.split(' ')[-2].rstrip('\n')
                    lastqueryindex = line.split(' ')[-1].rstrip('\n')
                    read_match = True
                    match_start = line.find('+') + 2
                    if new_query:
                        tmp = line.replace(" ", "")
                        query_start = tmp[3:tmp.find('+')]
                        new_query = False
                        #print "qstart"

                elif read_match:
                    #print(4)
                    #print line
                    match_seq += line[match_start:].rstrip('\n')
                    read_match = False

                elif line[:3] == "Tgt" and current_entry == target_lines[current_target_index]:
                    #print(5)
                    #print line
                    target_seq += line.split(' ')[-2].rstrip('\n')
                    lasttargetindex = line.split(' ')[-1].rstrip('\n')
                    if new_target:
                        tmp = line.replace(" ", "")
                        target_start = tmp[3:tmp.find('+')]
                        new_target = False

                elif line.find(" cols, ") > 0 and line.find(" ids ") > 0 and line.find(" gaps ") > 0:
                    if current_entry == target_lines[current_target_index]:
                        #print("writing...")
                        target_end = lasttargetindex
                        query_end = lastqueryindex
                        new_query = True
                        new_target = True
                        output.write(">%s\t%s\t%s\t%s\t%s\t%s\n" % (
                            query_name, target_name, query_start, query_end, target_start, target_end))
                        output.write("%s\n" % query_seq)
                        output.write("%s\n" % target_seq)
                        output.write("%s\n" % match_seq)
                        output.write("\n")
                        current_target_index = current_target_index + 1
                        #print "end %d" % current_entry
                    current_entry = current_entry + 1
                else:
                    pass

def screen(vsearch_aln_out_file, vsearch_caln_userout_file):
    # output_file = sys.argv[3]
    output_file = "%s.bad" % vsearch_aln_out_file

    frameshift_seqs, total_lines = flag_frame_shift_seqs(vsearch_caln_userout_file)
    frameshift_count = len(frameshift_seqs)
    print "Found frameshifts in %d / %d sequences" % (frameshift_count, total_lines)
    print "Gathering output data..."
    # print frameshift_seqs
    # vsearchLogFile = "test.log"
    # outputFile = "test.rslt"
    # bads = [i for i in range(10)]
    if len(frameshift_seqs) > 0:
        scrape_ulog(vsearch_aln_out_file, output_file, frameshift_seqs)