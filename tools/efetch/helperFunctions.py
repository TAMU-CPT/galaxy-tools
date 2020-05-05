import subprocess
import os

def cat_files(file_ext, output="cat_combination.txt",galaxy=False):
    """ concatenates files together """
    #cwd = os.getcwd()
    try:
        cmd = "cat *.{} > {}".format(file_ext,output)
        subprocess.run(cmd, shell=True)#, stdout=output)
    except subprocess.TimeoutExpired as err:
        print(err)
        
def awk_files(file_ext, output="awk_combination.txt",galaxy=False):
    try:
        cmd = "awk NF *.{} > {}".format(file_ext, output)
        subprocess.run(cmd, shell=True)
    except subprocess.TimeoutExpired as err:
        print(err)
def clean_files(input="combination.txt"):
    print("START HERE WITH REMOVING WHITE SPACE FROM ODD ADDITIONS")

if __name__ == "__main__":
    #cat_files("fasta")
    #clean_files()
    awk_files("DAT")