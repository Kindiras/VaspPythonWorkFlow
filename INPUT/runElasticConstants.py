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
    os.mkdir('nostrain')
    os.chdir('nostrain')
    incar.write_file('INCAR')
    kpoints.write_file('KPOINTS')
    structure.to(filename='POSCAR')
    pv.write_potcar(specs_file)
    pv.dumpFile('job_run_script', jobScript)
    subprocess.run('sbatch job_run_script',shell = True)
    os.chdir('../')
    
    delta_specs = specs_file['delta_list']
    change = np.linspace(delta_specs['begin'], delta_specs['end'], delta_specs['sample_points'])
    test_type_list,strain_matrix = pv.get_test_type_strain_delta_list()
    
    for test_type,strain in zip(test_type_list,strain_matrix):
        os.mkdir(str(test_type))
        os.chdir(str(test_type))
        for delta in change:
            os.mkdir(str(round(delta,4)))
            os.chdir(str(round(delta,4)))
            incar.write_file('INCAR')
            kpoints.write_file('KPOINTS')
            copystruc = structure.copy()
            copystruc.modify_lattice(strain(round(delta,4)))
            copystruc.to(filename='POSCAR')
            pv.write_potcar(specs_file)
            pv.dumpFile('job_run_script', jobScript)
            subprocess.run('sbatch job_run_script',shell = True)
            os.chdir('../')
        os.chdir('../')