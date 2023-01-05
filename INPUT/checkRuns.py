import os
import subprocess
import shutil
import glob



def readStdout(filename):
    with open(filename, "r") as file:
        last_line = file.readlines()[-1]
        last_line_words = [item for item in last_line.split()]
    return last_line_words    



def runVasp(N,paths,kdir,edir):
    if not os.path.isdir(paths+kdir+"/"+edir+"/run_"+str(N)):
        os.mkdir("run_"+str(N))

    odrive = paths+kdir+"/"+edir+"/"
    ndrive = paths+kdir+"/"+edir+"/run_"+str(N)
    output = ["POSCAR","CONTCAR", "POTCAR", "KPOINTS", "INCAR", "job_run_script"]
    files = glob.glob(odrive+'*')
    for f in files:
        outfiles = os.path.split(f)
        ovasp = (outfiles[1])
        if (ovasp in output):
            print("copying file",f)
            shutil.copy2(f,ndrive)
        else:
            print("moving file",f)
            shutil.move(f,ndrive)
    if(os.path.getsize("CONTCAR")>0):
        os.remove("POSCAR")
        shutil.move("CONTCAR","POSCAR")
    else:
        print("CONTCAR was not created from earlier run")       
    subprocess.call('sbatch job_run_script',shell=True) 



if __name__ == '__main__':
    N = 5
    paths = "/scratch/ikhatri/TaC_B1_SuperCells/3By3By3Tatom216/OPTM_ENCUTKPOINTS/"
    kdir_list = os.listdir(paths)
    for kdir in kdir_list:
        os.chdir(paths+kdir)
        print("checking on:",kdir)
        edir_list = os.listdir()
        for edir in edir_list:
            os.chdir(paths+kdir+"/"+edir)
            print("Working on:",edir)
            lastLineWords = readStdout("stdout")
            if("reached" in lastLineWords):
                print("It's already converged")
            else:    
                runVasp(N, paths, kdir, edir)
                print("Submitted in",edir)
            os.chdir("../")
        os.chdir("../")        
