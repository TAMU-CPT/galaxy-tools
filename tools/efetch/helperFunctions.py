import subprocess
import os

def cat_files(file_ext, direc=False, output="cat_combination.txt",galaxy=False):
    """ concatenates files together """
    #cwd = os.getcwd()
    try:
        if direc:
            os.chdir(direc)
            cmd = "cat *.{} > {}".format(file_ext,output)
        else:
            cmd = "cat *.{} > {}".format(file_ext,output)
        subprocess.run(cmd, shell=True)#, stdout=output)
    except subprocess.TimeoutExpired as err:
        print(err)


def awk_files(file_ext, output="awk_combination.txt", galaxy=False):
    try:
        if galaxy:
            cmd = "awk NF *.{} > {}Multi.{}".format(file_ext, output, file_ext)
        else:
            cmd = "awk NF *.{} > {}".format(file_ext, output)
        subprocess.run(cmd, shell=True)
    except subprocess.TimeoutExpired as err:
        print(err)


def pass_flag(input,flag="--output"):
    try:
        cmd = "{} {}".format(flag,input)
        subprocess.run(cmd,shell=True)
    except subprocess.TimeoutExpired as err:
        print(err)

def redirect(input):
    pass


if __name__ == "__main__":
    #cat_files("fasta")
    #clean_files()
    awk_files("DAT")