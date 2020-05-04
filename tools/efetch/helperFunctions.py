import subprocess
import os

def cat_files(file_ext):
    """ concatenates files together """
    cwd = os.getcwd()
    print(cwd)
    try:
        result = subprocess.run(["cat"] + [cwd+"/*."+file_ext] + ["> output.dat"])
    except subprocess.TimeoutExpired as err:
        result = err
    
    return result


if __name__ == "__main__":
    cat_files("fasta")