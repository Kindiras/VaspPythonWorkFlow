import argparse
import json
import yaml
import sys
import os
import numpy as np 
import re
import shutil
import subprocess
import matplotlib
matplotlib.use('Agg')
import datetime
import pymatgen as mg
from pymatgen.io import vasp
from scipy.optimize import curve_fit as cf
import matplotlib.pyplot as plt


# environmental variable or manual setting these lines required
VASP_EXEC = os.getenv('VASP_EXECC', 'mpirun -np 16 /home/ikhatri/local/vasp.std.hopper')
#VASP_EXEC = os.getenv('VASP_EXEC', 'OR-PATH-TO-YOUR-VASP-EXEC-default')
VASP_POTENTIALS_DIR = os.getenv('VASP_POTENTIALS_DIR', 'OR-PATH-TO-YOUR-VASP_POTENTIALS_DIR-default')
VASP_TEMPLATES_DIR = os.getenv('VASP_TEMPLATES_DIR', os.path.join(os.getcwd(), 'INPUT/TEMPLATES'))

def argpass():
    fileName = np.zeros(2,dtype = object)
    """
    This will pass the file name to the python.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('jobscript')
    args = parser.parse_args()
    fileName[0] = args.filename
    fileName[1] = args.jobscript
    return fileName

def fileload(filename):
    """

    Load a json or specs file, determined by the extension.

    """

    with open(filename, 'r') as f:
        if filename.endswith('.json'):
            file_dict = json.load(f)
        elif filename.endswith('.yaml'):
            file_dict = yaml.load(f,Loader=yaml.FullLoader)
    return file_dict
    

def filedump(dict_to_file, filename):
    """

    Dump a json or specs file, determined by the extension. Indentation of json
    and flow style of yaml is set.

    """

    with open(filename, 'w') as f:
        if filename.endswith('.json'):
            json.dump(dict_to_file, f, indent=4)
        elif filename.endswith('.yaml'):
            yaml.dump(dict_to_file, f, default_flow_style=False)



def get_structure(run_specs):
    """

    Get pymatgen.Structure. There are many ways to get a structure. They are all
    specified under 'poscar' tag in the specs file. There are two ways to get a
    structure.

    1. From an already made structure, either from a template POSCAR, or a
    Materials Project database entry. (Recommended)

    If 'template' is present in 'poscar', you can either set the
    VASP_TEMPLATES_DIR environmental variable, or just leave it to the default
    'INPUT/TEMPLATES'. After that, this specified POSCAR-type template file path
    will be used to obtain the structure from the VASP_TEMPLATES_DIR, and a
    structure is returned. If you set 'rel_to_run_dir' to True, 'template'
    refers to the file relative to the 'run_dir'.

    If 'material_id' is present in 'poscar', MAPI_KEY environmental variable
    needs to be set according to the Materials Project (materialsproject.org).

    An optional 'get_structure' can be set to one of

        ['primitive', 'reduced', 'sorted', 'conventional_standard',
         'primitive_standard', 'refined']

    For 'primitive', 'conventional_standard' and so on, an additional tag 'prec'
    controls the tolerence/symmetry finding precision threshold.

    The 'primitive', 'reduced' and 'sorted' are methods of the object
    pmg.Structure, while the rest are methods of
    pmg.symmetry.analyzer.SpacegroupAnalyzer. Please refer to the pymatgen docs
    to see what they are exactly. Be careful about whether the resultent
    structure is what you want. When in doubt, always manually run the pymatgen
    commands and see the outcome.

    By default, the code uses the element types written in the structure to
    generate the POTCAR by maintaining the existence of 'elem_types'. However,
    if you set 'repl_elems' with a dict, like {N: C, Ti: Ti_sv}, the elements in
    the structure will be accordingly replaced and the necessary mechanism is in
    place to make sure POTCAR is to be generated with the flavored potentials.

    Setting 'elem_types' in the specs file as a list of potentials (can have
    flavors like Ti_sv) in the same sequence as in the structure also works, but
    one has to be careful to match the sequence correctly. Not recommended.

    An optional 'volume' can be set to scale the structure of the template.

    2. From manual description. (Cumbersome)

    The manual generation from spacegroup is done by specifying

    'spacegroup' (international number or symbol)

    'cryst_sys' (one of the seven)

    'lattice_params' ('a', 'b', 'c', 'alpha', 'beta', 'gamma', some of which
    'can be omitted because of a more symmetric crystallographic system)

    'elem_types' (the elements in the structure, which can be flavored
    'potentials, e.g. Ti_sv)

    'atoms_multitude' (multitude of atoms of the same element type in a list,
    'the sequence following 'elem_types'. Only symmetrically distinct species and
    'coords should be provided, according to the Wychoff positions)

    'atoms_direct_coords' (direct locations, relative to the lattice vectors
    'of the symmetrically distinct atoms. There should be the same number of
    'them as the sum of atoms_multitude)

    """

    is_template = None
    is_material_id = None
    poscar_specs = run_specs['poscar']

    if 'template' in poscar_specs:
        is_template = True
        if 'rel_to_run_dir' in poscar_specs and poscar_specs['rel_to_run_dir']:
            poscar = mg.io.vasp.Poscar.from_file(poscar_specs['template'])
        else:
            poscar = mg.io.vasp.Poscar.from_file(os.path.join(VASP_TEMPLATES_DIR, poscar_specs['template']))
        structure = poscar.structure
    elif 'material_id' in poscar_specs:
        is_material_id = True
        m = mg.MPRester()
        structure = m.get_structure_by_material_id(poscar_specs['material_id'])
    if is_template or is_material_id:
        if 'get_structure' in poscar_specs:
            prec = poscar_specs['prec'] if 'prec' in poscar_specs else 1e-3
            sga = mg.symmetry.analyzer.SpacegroupAnalyzer(structure, symprec=prec)
            if poscar_specs['get_structure'] == 'sorted':
                structure = structure.get_sorted_structure()
            if poscar_specs['get_structure'] == 'reduced':
                structure = structure.get_reduced_structure()
            elif poscar_specs['get_structure'] == 'primitive':
                structure = structure.get_primitive_structure(prec)
            elif poscar_specs['get_structure'] == 'primitive_standard':
                structure = sga.get_primitive_standard_structure()
            elif poscar_specs['get_structure'] == 'conventional_standard':
                structure = sga.get_conventional_standard_structure()
            elif poscar_specs['get_structure'] == 'refined':
                structure = sga.get_refined_structure()

        symbol_set = list(structure.symbol_set)
        
        if 'repl_elems' in run_specs:
            for idx, symbol in enumerate(symbol_set):
                if symbol in run_specs['repl_elems']:
                    symbol_set[idx] = run_specs['repl_elems'][symbol]
                    #print(list(run_specs['repl_elems'].values()))
            repl_elems_struct = {key: re.sub(r'_.*', '', value) for key, value in run_specs['repl_elems'].items()}
            structure.replace_species(repl_elems_struct)
        run_specs['elem_types'] = list(run_specs['repl_elems'].values())

        if 'volume' in poscar_specs:
            structure.scale_lattice(poscar_specs['volume'])
        return structure
    else:
        cryst_sys = poscar_specs['cryst_sys']
        lattice_params = poscar_specs['lattice_params']
        if cryst_sys == 'cubic':
            lattice = mg.Lattice.cubic(lattice_params['a'])
        elif cryst_sys == 'hexagonal':
            lattice = mg.Lattice.hexagonal(lattice_params['a'], lattice_params['alpha'])
        elif cryst_sys == 'tetragonal':
            lattice = mg.Lattice.tetragonal(lattice_params['a'], lattice_params['c'])
        elif cryst_sys == 'orthorhombic':
            lattice = mg.Lattice.orthorhombic(lattice_params['a'], lattice_params['b'], lattice_params['c'])
        elif cryst_sys == 'rhombohedral':
            lattice = mg.Lattice.rhombohedral(lattice_params['a'], lattice_params['alpha'])
        elif cryst_sys == 'monoclinic':
            lattice = mg.Lattice.orthorhombic(lattice_params['a'], lattice_params['b'], lattice_params['c'],
                lattice_params['beta'])
        else:
            lattice = mg.Lattice.orthorhombic(lattice_params['a'], lattice_params['b'], lattice_params['c'],
                lattice_params['alpha'], lattice_params['beta'], lattice_params['gamma'])

        elem_types_struct = [re.sub(r'_.*', '', i) for i in poscar_specs['elem_types']]
        elem_types_struct_multi = []
        for i, elem in enumerate(elem_types_struct):
            elem_types_struct_multi.extend([elem] * poscar_specs['atoms_multitude'][i])

        structure = mg.Structure.from_spacegroup(poscar_specs['spacegroup'], lattice,
                elem_types_struct_multi, poscar_specs['atoms_direct_coords'])
        run_specs['elem_types'] = poscar_specs['elem_types']
        return structure




def read_incar(run_specs):
    """

    Read contents of 'incar' from the specs file. If 'incar' does not exist,
    return an empty Incar object, still functional.

    """

    incar = mg.io.vasp.Incar()
    if 'incar' in run_specs and run_specs['incar']:
        incar.update(run_specs['incar'])
    return incar




def read_kpoints(run_specs, structure=None):
    """

    Read contents of 'kpoints' from the specs file. If 'kpoints' does not
    exist, return an automatic (A) mesh with 5 subdivisions along each
    reciprocal vector.

    If 'kpoints' exists, it follows the sequence below.

    If 'density' exists in 'kpoints', automatic density mesh is returned, with
    density being KPPRA. In this case you need to provide the structure as an
    argument. 'force_gamma' can be specified to True.

    Otherwise, if 'divisions', specified as a list, exists in 'kpoints', MP mesh or Gamma-centered
    mesh is returned, according to 'mode' that starts with 'M' or 'G'.

    """

    kpoints = mg.io.vasp.Kpoints.automatic(5)
    if 'kpoints' in run_specs and run_specs['kpoints']:
        kpoints_spec = run_specs['kpoints']
        if 'density'in kpoints_spec:
            if 'force_gamma'in kpoints_spec:
                force_gamma = kpoints_spec['force_gamma']
            else:
                force_gamma = False
            kpoints = mg.io.vasp.Kpoints.automatic_density(structure, kpoints_spec['density'],
                force_gamma=force_gamma)
        elif 'divisions' in kpoints_spec:
            if kpoints_spec['mode'].upper().startswith('M'):
                kpoints = mg.io.vasp.Kpoints.monkhorst_automatic(kpoints_spec['divisions'])
            elif kpoints_spec['mode'].upper().startswith('G'):
                kpoints = mg.io.vasp.Kpoints.gamma_automatic(kpoints_spec['divisions'])
    return kpoints

def get_run_dir(run_specs):
    """

    Get the directory where the routine takes place.

    If 'run_dir' is in the specs file, use that; otherwise name it 'vasp_test'.

    """

    dirname = run_specs['run_dir'] if 'run_dir' in run_specs else 'vasp_test'
    return dirname



def chdir(dirname):
    """

    Enter a path. If it does not exist, create one recursively and enter.

    """

    os.makedirs(dirname, exist_ok=True)
    os.chdir(dirname)    



def write_potcar(run_specs):
    """

    Write POTCAR. It gets the POTCAR types from 'elem_types' in the specs file.
    You need to set the VASP_POTENTIALS_DIR environmental variable, or edit
    the head of this file to be able to use it.

    """

    potential_base = os.path.join(VASP_POTENTIALS_DIR, run_specs['pot_type'])
    with open('POTCAR', 'wb') as outfile:
        for filename in [os.path.join(potential_base, e, 'POTCAR') for e in run_specs['elem_types']]:
            with open(filename, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)


def fileLoad(filename):
    with open(filename,"rt") as f:
            data = [line.rstrip() for line in f]
        #data = [[i for i in line.split()] for line in f.readlines()]
    return data



def dumpFile(filename,inData):
    with open(filename,'w') as fp:
        for lines in inData:
            fp.write(str(lines))
            fp.write("\n")




"""def run_vasp(filename):
    time_format = ' "\n----------\nreal     %E" '
    time = '/usr/bin/time -f ' + time_format
    print("time:",time)
    arg = 'sbatch'+' '+ filename
    print(subprocess.check_call(time + VASP_EXEC + ' 2>&1 | tee -a stdout', shell=True))
    p1=subprocess.call(arg, stdout=subprocess.PIPE, shell=True)
    print( "hbreak",' "\n' + '=' * 100 + '\n" ')  
    #process = subprocess.call(arg,shell=True)
    print(p1)   """         





def get_test_type_strain_delta_list():
    
    strain_matrix = [] # 3x3 strain matrix D => a = aD
    test_type_list = ["c11+2c12", "c11-c12", "c44"]

    strain_matrix.append(lambda delta: np.array([[1 + delta, 0, 0],
                                                   [0, 1 + delta, 0],
                                                   [0, 0, 1 + delta]]))

    strain_matrix.append(lambda delta: np.array([[1 + delta, 0, 0],
                                                   [0, 1 - delta, 0],
                                                   [0, 0, 1 + delta ** 2 / (1 - delta ** 2)]]))

    strain_matrix.append(lambda delta: np.array([[1, delta/2, 0],
                                                   [delta/2, 1, 0],
                                                   [0, 0, 1 + delta ** 2 / (4 - delta ** 2)]]))

    return test_type_list, strain_matrix




def curve_fit(fit_eqn, X, Y, p0=None, sigma=None, absolute_sigma=False, plot=False, on_figs=None, **kw):
    """
    Fit the two equal length arrays to a certain function Y = f(X).

    Parameters
    ----------
    X: array
    Y: array
    plot: bool
        plot figures or not. Default to False.
    on_figs: list/int
        the current figure numbers to plot to, default to new figures

    Returns
    -------
    a dict, containing
        'params': fitting parameters
        'r_squared': value for evaluating error
        'fitted_data': a dict that has 2D array of fitted data
            easily to Pandas DataFrame by pd.DataFrame(**returned_dict['fitted_data'])
        'ax': the axes reference, if return_refs == True
    """
    X = np.array(X)
    Y = np.array(Y)
    popt, pcov = cf(fit_eqn, X, Y, p0, sigma, absolute_sigma, **kw)
    X_fit = np.linspace(sorted(X)[0], sorted(X)[-1], 1000)
    Y_fit = fit_eqn(X_fit, *popt)
    Y_fit_eqlen = fit_eqn(X, *popt)
    r_squared = get_r_squared(Y, Y_fit_eqlen)
    data = np.column_stack((X_fit, Y_fit))
    col_names = ['X_fit', 'Y_fit']

    return_dict = {'params': popt, 'r_squared': r_squared,
        'fitted_data': {'columns': col_names, 'data': data}}

    if plot:
        initiate_figs(on_figs)
        plt.plot(X, Y, 'o')
        plt.plot(X_fit, Y_fit, '-')
        axes = {'ax': plt.gca()}
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.tight_layout()
        return_dict.update(axes)

    return return_dict



def initiate_figs(on_figs):
    if on_figs is None:
        plt.figure()
    elif isinstance(on_figs, int):
        plt.figure(on_figs)
    else:
        plt.figure(on_figs.pop(0))




def central_poly(X, a, b, c):
    return b * X**3 + a * X**2 + c



def get_r_squared(Y, Y_fit_eqlen):
    ss_resid = np.sum((Y_fit_eqlen - Y) ** 2)
    Y_avg = np.sum(Y) / len(Y)
    ss_total = np.sum((Y - Y_avg) ** 2)
    r_squared = 1 - ss_resid / ss_total
    return r_squared



def solve(combined_econst_array):
    """

    Solves for the elastic constants from the matrix and coeffs.


    Used in combination with the function get_test_type_strain_delta_list. The
    argument combined_econst_array is an array of the second order polynomial
    coefficients coming from the runs with various test_type, in the same order.

    """
    econsts_str = ["C11", "C12", "C44"]
    coeff_matrix = np.array([[3/2., 3, 0],
                                 [1, -1, 0],
                                 [0, 0, 1/2.]])

    solved = np.linalg.solve(coeff_matrix, combined_econst_array)
    return dict(zip(econsts_str, solved))



def run_vasp():
    """

    Run VASP, time it, print out to screen, and append it to a file named
    'stdout'. You need to set the VASP_EXEC environmental variable, or edit
    the head of this file to be able to use it.

    """

    time_format = ' "\n----------\nreal     %E" '
    time = '/usr/bin/time -f ' + time_format
    subprocess.check_call(time + VASP_EXEC + ' 2>&1 | tee -a stdout', shell=True)
    hbreak = ' "\n' + '=' * 100 + '\n" '
    subprocess.call('echo -e ' + hbreak + ' | tee -a stdout', shell=True)

    # print(subprocess.check_output(time + VASP_EXEC, shell=True))
    # time = timeit.Timer('print(subprocess.check_output("{}", shell=True).decode())'.format(VASP_EXEC), 'import subprocess').timeit(1)
    # print('\n----------\nreal     ' + str(datetime.timedelta(seconds=int(time))))
    # print('\n' + '=' * 100 + '\n')

