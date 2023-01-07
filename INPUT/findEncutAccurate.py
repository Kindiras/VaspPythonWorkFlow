import matplotlib.pyplot as plt
import pymatgen as mg
from  pymatgen.io import vasp
import os
import numpy as np

paths = "/scratch/ikhatri/HfC/ConverntionalCell"
os.chdir(paths)
#os.mkdir("RESULTS")
def sort_numeric_strings(nums_str):
    result = sorted(nums_str, key=lambda el: int(el))
    return result




Energy = []
Encut = []
Eref = 0
kdirs = os.listdir(paths+"/OPTM_ENCUT")
kdir =sort_numeric_strings(kdirs)

for idx,dir in enumerate(kdir):
    print("working on:",dir)
    os.chdir(paths+"/OPTM_ENCUT/"+dir)
    edirs = os.listdir()
    edir = sort_numeric_strings(edirs)
    for edr in edir:
        os.chdir(paths+"/OPTM_ENCUT/"+dir+"/"+edr)
        oszicar = mg.io.vasp.Oszicar('OSZICAR')
        poscar =mg.io.vasp.Poscar.from_file("POSCAR")
        st = poscar.structure
        Energy.append(round(oszicar.final_energy,8))
        Encut.append(int(edr))
        os.chdir('../')

plt.plot(Encut, Energy, '*')
plt.xlabel('ENCUT (eV)')
plt.ylabel('Energy (eV)')
plt.tight_layout()
plt.savefig(paths+"/RESULTS/energy-encut.pdf")
plt.close()
np.savetxt(paths+"/RESULTS/energy-encut.txt", np.column_stack((Encut, Energy)), '%12.4f', header='encut energy')

    