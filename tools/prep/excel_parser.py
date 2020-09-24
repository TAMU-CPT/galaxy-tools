import re

import pandas as pd
#from pyxlsb import open_workbook


class ExcelParsing:
    
    def __init__(self, file):
        self.file = file
    

    @staticmethod
    def parse_excel(file):
        """ returns a dataframe object, checks binary """
        
        #try:
        return pd.read_excel(file)
        #except:
        #    df=[]
        #    with open_workbook(file) as wb:
        #        with wb.get_sheet(1) as sheet:
        #            for row in sheet.rows():
        #                df.append([item.v for item in row])
        #    return pd.DataFrame(df[1:], columns=df[0])
    
    def chop_frame(self, cols=None, **kwargs):
        """ subset based on columns """

        frame = self.parse_excel(self.file, **kwargs)
        return frame.loc[:,cols]


if __name__ == "__main__":
    
    file = 'test-data/Refseq_coliphages_20200908.xlsb'
    data = ExcelParsing(file).chop_frame(cols=['Genome','Accession'])
    names = list(data['Genome'])
    spliced_names = []
    for name in names:
        s = name.split(' ')
        if re.search(('phage|virus|coli'), s[1]):
            r = s[2:]
        else:
            r = s[1:]
        spliced_names.append(r)
    print(spliced_names)
    print(len(names))
    print(len(spliced_names))
        #if name.count(' ') == 2:
        #    print(f'{name} has 2')
        #else:
        #    print(f'{name} ha XX')
    #print(list(data['Genome'].str.split(' ').str[-1]))