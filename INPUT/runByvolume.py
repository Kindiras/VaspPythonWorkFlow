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
    volume_specs = specs_file['volume_test']
    change = np.linspace(volume_specs['begin'],volume_specs['end'],volume_specs['sample_points'])
    for volumes in change:
        os.mkdir(str(round(volumes,4)))
        os.chdir(str(round(volumes,4)))
        structure.scale_lattice(volumes)
        incar.write_file('INCAR')
        kpoints.write_file('KPOINTS')
        structure.to(filename='POSCAR')
        pv.write_potcar(specs_file)
        pv.dumpFile('job_run_script',jobScript)
        subprocess.run('sbatch job_run_script',shell=True)
        os.chdir('../')