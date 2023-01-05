#from turtle import color
import matplotlib.pyplot as plt
import pymatgen as mg
from  pymatgen.io import vasp
import os
import numpy as np
import pandas as pd

paths = "//scratch/ikhatri/TaC_B1_SuperCells/3By3By3Tatom216"
os.chdir(paths)
os.mkdir("RESULTS")
def sort_numeric_strings(nums_str):
    result = sorted(nums_str, key=lambda el: int(el))
    return result

rows = int(input("how many kpoints ?"))
cols = int(input("how many Encuts ?"))

Energy = [[0]*cols for i in range(rows)]
kpointNum = np.zeros(rows)
encutNum = np.zeros(cols)
kPos = 0

kdir = os.listdir(paths+"/OPTM_ENCUTKPOINTS")
kpdirs = sort_numeric_strings(kdir)
for kdr in kpdirs:
    kpointNum[kPos] = int(kdr)
    print("Working on Kpoints of ",kdr)
    print("***************************************")
    os.chdir(paths+"/OPTM_ENCUTKPOINTS/"+kdr)
    enPos = 0
    edir = os.listdir()
    endirs = sort_numeric_strings(edir)
    for endr in endirs:
        print("Working on Encut ",endr)
        os.chdir(paths+"/OPTM_ENCUTKPOINTS/"+kdr+"/"+endr)
        oszicar = mg.io.vasp.Oszicar('OSZICAR')
        poscar =mg.io.vasp.Poscar.from_file("POSCAR")
        st = poscar.structure
        Energy[kPos][enPos] = oszicar.final_energy/st.num_sites
        encutNum[enPos] = int(endr)
        enPos = enPos + 1
        os.chdir('../')
    kPos = kPos + 1    
    os.chdir('../')   

colors = ['black','gray','gold','green','red','chocolate']
#colors = ['black','gray','gainsboro','maroon','red','chocolate','peru','gold','lawngreen','green','cyan','blue']
#plt.plot(EK,Energy)
#plt.savefig("/scratch/ikhatri/DFT_TaC/TaC_ZB/RESULTS/EnergyEncutKpoint.pdf")
#plt.close()
for i,array in enumerate(Energy):
    plt.plot(encutNum,array,color=colors[i],marker = "o",label = f"KPOINTS {int(kpointNum[i])}")
    #plt.plot(encutNum,array,color=np.random.rand(3,),marker = "o",label = f"KPOINTS {int(kpointNum[i])}")
plt.legend(loc = "lower center",bbox_to_anchor=(0.5,0.40))    
plt.savefig(paths+"/RESULTS/EncutKpoitTest.pdf")
plt.close()
Energies = np.array(Energy).T
DF = pd.DataFrame(Energies,columns=[i for i in kpointNum],index=[j for j in encutNum])
DF.to_csv(paths+"/RESULTS/EncutKpoitTest.csv")
