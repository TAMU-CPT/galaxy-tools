from time import sleep
from Bio import Entrez
from Bio import SeqIO


class CPTLink:

    def __init__(self, email, db, acc, dbfrom=None):
        self.email = email
        self.db = db
        self.dbfrom = dbfrom
        self.acc = acc
        Entrez.email = self.email

    def __repr__(self):
        return f"< email: {self.email} | accession: {self.acc} | database: {self.db}>"
    
    @staticmethod
    def grab_gacc_with_gid(db, gid):
        """
        Grabs Genome accession with NCBI gid
        """
        esum = Entrez.read(Entrez.esummary(db=db,id=gid))
        return esum[0]["AccessionVersion"]

    @staticmethod
    def get_tax_with_gacc(dbfrom,acc,tax='taxonomy'):
        """
        Uses genome accession and retrieves tax id
        db is defaulted to taxonomy
        """
        tax_link = Entrez.read(Entrez.elink(db=tax,dbfrom=dbfrom,id=acc))
        tax_id = tax_link[0]['LinkSetDb'][0]['Link'][0]['Id']
        
        ui_list = Entrez.read(Entrez.efetch(db='taxonomy', id=tax_id, ret_type='uilist', ret_mode='xml'))
        org = ui_list[0]['ScientificName']

        return tax_id,org

    def map_accessions(self,wp_all=False):
        """
        Maps accessions from an input database to another database. Currently works best if 
        mapping from the 'protein' database to the 'nuccore' database. 
        """
        try:
            record = Entrez.read(Entrez.elink(db=self.db, dbfrom=self.dbfrom, id=self.acc))
        except RuntimeError:
            sleep(15)
            record = Entrez.read(Entrez.elink(db=self.db, dbfrom=self.dbfrom, id=self.acc))
        sleep(1.05) #  don't overload NCBI.
        if self.acc[0:2] == "WP": #  do some different stuff with WPs
            WP_amt = len(record[0]['LinkSetDb'][0]['Link'])
            if wp_all:
                list_of_WPs = []
                for each_mapping_id in record[0]['LinkSetDb'][0]['Link']:
                    acc = self.grab_gacc_with_gid(self.db,each_mapping_id['Id'])
                    tax_id, org = self.get_tax_with_gacc(self.db,acc)
                    list_of_WPs.append([acc,tax_id,org])

                return self.acc, list_of_WPs, WP_amt
            else:
                gid = record[0]['LinkSetDb'][0]['Link'][0]['Id']
                acc = self.grab_gacc_with_gid(self.db,gid)
                tax_id, org = self.get_tax_with_gacc(self.db,acc)

                return self.acc, [[acc,tax_id,org]], WP_amt
        else:
            raise "Not a WP accession"


if __name__ == "__main__":
    pass
