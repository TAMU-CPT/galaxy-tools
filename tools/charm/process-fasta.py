import sys


if __name__ == '__main__':
    seen = False
    with open(sys.argv[1], 'r') as handle:
        for line in handle:
            if line.startswith('>'):
                if seen:
                    sys.exit(0)
                seen = True
                print '>' + sys.argv[2]
            else:
                print line,
