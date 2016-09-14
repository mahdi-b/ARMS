import sys
def check(inputFile):
    for line in open(inputFile, 'r' ):
        data = line.split(' ')
        total = data[0].split("_")[-1]
        print "_".join(data[0].split("_")[:-1])
        sum = 0
        for item in data[1:]:
            sum += int(item.split("_")[-1])
        if total == sum:
            print "OK"
        else:
            print "Found %d, expected %d" % (int(sum), int(total))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: name_file"
        exit()
    check(*sys.argv[1:2])