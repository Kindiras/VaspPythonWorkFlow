from genericpath import isfile
import subprocess
import pyVaspWorkFlow as pv
import os
import numpy as np
import matplotlib.pyplot as plt
import pymatgen as mg
from pymatgen.io import vasp


if __name__ == '__main__':

    file_name = pv.argpass()
    specs_file = pv.fileload(file_name[0])
    jobScript = pv.fileLoad(file_name[1])
    pv.chdir(pv.get_run_dir(specs_file))

    structure = pv.get_structure(specs_file)
    incar = pv.read_incar(specs_file)
    kpoints = pv.read_kpoints(specs_file,structure)

    encut_specs = specs_file['encut_test']
    change = np.array(range(encut_specs['begin'], encut_specs['end'], encut_specs['step']))

    kpoints_specs = specs_file['kpoints_test']
    if 'force_gamma' in kpoints_specs:
        force_gamma = kpoints_specs['force_gamma']
    else:
        force_gamma = False    
    density_change = np.array(kpoints_specs['density_change'])

    for kp in density_change:
        a = int(((kp*structure.lattice.a*structure.lattice.b*structure.lattice.c)/structure.num_sites)**(1/3)/structure.lattice.a)
        c = int(((kp*structure.lattice.a*structure.lattice.b*structure.lattice.c)/structure.num_sites)**(1/3)/structure.lattice.c)
       
        if os.path.isdir(str(a)+str(c)):
            os.chdir(str(a)+str(c))
        else:
            os.mkdir(str(a)+str(c))
            os.chdir(str(a)+str(c))
        #root_cube = lambda x: x**(1./3.)#This is for cubic structure B1,B2 or B3
        #dirName = int(root_cube(kp/structure.num_sites))
        #os.mkdir(str(dirName))
        for encut in (change):
            os.mkdir(str(encut))
            os.chdir(str(encut))
            incar['ENCUT'] = encut
            incar.write_file('INCAR')
            kpoints = vasp.Kpoints.automatic_density(structure, kp, force_gamma=force_gamma)
            kpoints.write_file('KPOINTS')
            structure.to(filename='POSCAR')
            pv.write_potcar(specs_file)
            pv.dumpFile('job_run_script', jobScript)
            subprocess.run('sbatch job_run_script',shell = True)
            os.chdir('../')
        os.chdir('../')    