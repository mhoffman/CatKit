import importlib
import json
import copy
import sys
import ase

from . import fwio

try:
    import espresso
except(ImportError, KeyError):
    espresso = None

def str_to_class(full_classname):
    mod_name, classname =  full_classname.rsplit('.', 1)
    return getattr(importlib.import_module(mod_name),
                  classname)


def get_potential_energy(
        calculator_name,
        in_file='input.traj',
        out_file='pw.pwo'):
    """Performs a ASE get_potential_energy() call with a compatible ase calculator.
    Keywords are defined inside the atoms object information.

    This can be a singlepoint calculation or a full relaxation depending
    on the keywords.

    Parameters
    ----------
    calculator : str
        String representation of a calculator import
    in_file : str
        Name of the input file to load from the local directory.
    out_file : str
        Name of the output file to read the completed trajectory form.
    """
    atoms = ase.io.read(in_file)

    # Setting up the calculator
    calculator = str_to_class(calculator_name)
    if calculator_name == 'espresso.Espresso':
        calc = calculator(atoms, **atoms.info['calculator_parameters'])

    else:
        calc = calculator(**atoms.info)
        atoms.set_calculator(calc)

    # Perform the calculation and write trajectory from log.
    atoms.get_potential_energy()

    if calculator_name == 'espresso.Espresso':
        # Patch for reading magmom of trajectory
        images = espresso.io.read(out_file, ':')
    elif calculator_name == 'ase.calculators.emt.EMT':
        images = copy.deepcopy(atoms)
    else:
        images = ase.io.read(out_file, ':')

    # Moneky patch for constraints
    for image in images:
        if hasattr(image, 'constraints'):
            image.constraints = atoms.constraints
        if hasattr(atoms, 'pbc') and hasattr(image, '_pbc'):
            image._pbc = atoms.pbc

    return fwio.atoms_to_encode(images)
