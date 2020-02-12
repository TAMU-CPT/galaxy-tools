############ Takes an input file, and converts the file to a list to be parsed

##### Should have a check and see if the file is "large"

class Synonyms:

    """
    This class reads a flat file as an object and parses it, based on input delimiters.
    """
    
    def __init__(self, file, delims='\n', comments='#'):
        self.file = file
        self.delims = delims
        self.comments = comments
        self.read = open(file, 'r')

    def parse_it(self):
        """
        iterates and prints the contents
        """
        with self.read as f:
            for item in f:
                print(item)
    
    def store_it(self):
        """
        stores items in a list
        """
        l = []
        if self.delims == '\n':
            with self.read as f:
                for line in f.readlines():
                    line = item.split()[0]
                    l.append(line)
            return l
        else:
            print("I can't parse this yet...")

if __name__ == "__main__":
    filename = 'tools/proximity/test-data/iterate.txt'
    p = Synonyms(filename,delims='\n')
    #p.parse_it()
    print('=====================================')
    l = p.store_it()
    print(l)
    for i in l:
        print(i)

    