from netrc import netrc
import importlib
import json
import pprint
try:  # Python 2/3 compatible import
    import io as StringIO
except ImportError:
    import StringIO

import ase.db.row
import fireworks
import numpy as np

from . import fwio
from . import fwase
import catkit.gen.surface
import catkit.gen.adsorption
import catkit.build


def open_string(string, format='cif'):
    with StringIO.StringIO() as temp_file:
        temp_file.write(string)
        temp_file.seek(0)
        atoms = ase.io.read(
            temp_file,
            format=format,
        )

    return atoms


def str_to_class(full_classname):
    mod_name, classname = full_classname.rsplit('.', 1)
    return getattr(importlib.import_module(mod_name),
                   classname)


def setup_calculator(atoms, calculator_name):
    # Setting up the calculator
    calculator = str_to_class(calculator_name)
    if calculator_name == 'espresso.Espresso':
        calc = calculator(atoms, **atoms.info.get('calculator_parameters', {}))

    else:
        calc = calculator(**atoms.info.get('calculator_parameters', {}))
    atoms.set_calculator(calc)


def calculate_bulk(
        calculator,
        bulk_params,
):
    print("============================")
    print("IN CALCULATE BULK")
    print(calculator)
    print(bulk_params)

    atoms = open_string(bulk_params['wyckoff']['cif'])
    setup_calculator(atoms, calculator)
    energy = atoms.get_potential_energy()

    print(energy)
    print("============================")
    return fwio.atoms_to_encode(atoms)


def calculate_empty_slab(
        calculator,
        slab_params,
        bulk_trajectory,
):

    print("============================")
    print("IN CALCULATE EMPTY SLAB")
    # print(calculator)
    # print(slab_params)
    # print(bulk_trajectory)
    # create empty slab

    bulk = fwio.encode_to_atoms(bulk_trajectory)[-1]
    slab = catkit.gen.surface.SlabGenerator(
        bulk,
        miller_index=(
            slab_params['millerX'],
            slab_params['millerY'],
            slab_params['millerZ'],
        ),
        layers=slab_params['layers'],
        fixed=slab_params['fixed'],
        vacuum=slab_params['vacuum'],
        standardize_bulk=True,
    ).get_slab()
    setup_calculator(slab, calculator)
    slab.get_potential_energy()

    print("============================")
    return fwio.atoms_to_encode(slab)


def calculate_adsorbates(
    calculator,
    slab_params,
    adsorbate_params,
    bulk_trajectory
):
    print("============================")
    print("IN CALCULATE ADSORBATE SLAB")
    # print(calculator)
    # print(slab_params)
    print(adsorbate_params)
    print(adsorbate_params.keys())
    # print(bulk_trajectory)
    bulk = fwio.encode_to_atoms(bulk_trajectory)[-1]
    slab = catkit.gen.surface.SlabGenerator(
        bulk,
        miller_index=(
            slab_params['millerX'],
            slab_params['millerY'],
            slab_params['millerZ'],
        ),
        layers=slab_params['layers'],
        fixed=slab_params['fixed'],
        vacuum=slab_params['vacuum'],
        standardize_bulk=True,
    ).get_slab()

    adsorbate = catkit.build.molecule(
        adsorbate_params['adsorbate']
    )[0]

    builder = catkit.gen.adsorption.Builder(
        slab,
    )

    slabs = builder.add_adsorbate(
        adsorbate,
        index=-1,
        bonds=[0],
        # TODO: works for atomic adsorbates
        # should be generalized
    )

    for slab in slabs:
        setup_calculator(slab, calculator)
        slab.get_potential_energy()

    print("============================")
    return str([
        fwio.atoms_to_encode(slab)
        for slab in slabs
    ])


def calculate_molecules(
        calculator,
        adsorbate_params,
):
    print("============================")
    print("IN CALCULATE MOLECULES")
    # print(calculator)
    # print(adsorbate_params)

    print("============================")

    TODO CONTINUE HERE
    return []


def collect_energies(
    input_data,
    bulk_trajectory,
    slab_trajectory,
    adsorbate_trajectories,
    molecule_trajectories,
):
    print("============================")
    print("IN COLLECTING ENERGIES")
    # print(input_data)
    # print(bulk_trajectory)
    # print(slab_trajectory)
    # print(adsorbate_trajectories)
    # print(molecule_trajectories)

    with open('collected_energies.json', 'w') as outfile:
        outfile.write(json.dumps({
            'input_data': input_data,
            'bulk_trajectory': bulk_trajectory,
            'slab_trajectory': slab_trajectory,
            'adsorbate_trajectories': adsorbate_trajectories,
            'molecule_trajectories': molecule_trajectories,
        }))
        # outfile.write(json.dumps(kwargs))

    print("============================")
    return None


class Cohesive():
    """Simple submission script helper for CatFlow.
    """

    def __init__(
            self,
            host,
            username=None,
            name=None,
            password=None,
            calculator='espresso.Espresso'
    ):
        """Initialize a fireworks instance."""
        if username is None or name is None or password is None:
            username, name, password = netrc().authenticators(host)
        if username is None:
            raise ValueError('username, name, and password required.')

        launchpad = fireworks.LaunchPad(
            host=host,
            name=name,
            username=username,
            password=password)

        self.launchpad = launchpad
        self.calculator = calculator

    def submit_adsorption(
            self,
            image,
            input_data,
            calculation_name,
            parameters=None,
            spec=None):
        """
        Run a relaxation of a given DB entry or atoms object. If a
        database object is used, the calculation will automatically store
        the keys and data for later retrieval.

        The entries uuid will also be stored and `data.calculator_parameters`
        will be used as the calculation parameters.

        Parameter:
        ----------
        images : Atoms object | AtomsRow object
            ASE database entry or atoms object to relax.
        calculation_name : str
            Name of the fireworks calculation to be used.
        parameters : dict
            Calculation parameters to use. Will be pulled from
            a database entry `data.calculator_parameters`.
        spec : dict
            Additional fireworks specifications to pass to the database.
        """
        keys, data = {}, {}
        if isinstance(image, ase.db.row.AtomsRow):
            atoms = image.toatoms()
            keys.update(image.key_value_pairs)
            keys.update({'uuid': image.unique_id})
            data.update(image.data)
        else:
            atoms = image

        if parameters is not None:
            atoms.info = parameters
        elif data.get('calculator_parameters'):
            atoms.info = data.get('calculator_parameters')
            del data['calculator_parameters']
        elif atoms.info:
            pass
        # else:
            #raise ValueError('Calculation parameters missing.')

        for k, v in data.items():
            if isinstance(v, np.ndarray):
                fwio.array_to_list(v)
                data[k] = v

        # pprint.pprint(input_data)
        bulk_cif = input_data.get('bulkParams', {}) \
                             .get('wyckoff', {}) \
                             .get('cif', '')

        bulk = open_string(bulk_cif)
        bulk_encoding = fwio.atoms_to_encode(bulk)

        tasks = []

        encoding = fwio.atoms_to_encode(atoms)

        # Relax Bulk
        tasks.append(fireworks.PyTask(
            func='catkit.flow.cohesive.calculate_bulk',
            args=[
                self.calculator,
                input_data['bulkParams'],
            ],
            outputs=['bulk_trajectory'],
            auto_kwargs=True,
        ))

        # Cut Empty Slab
        tasks.append(fireworks.PyTask(
            func='catkit.flow.cohesive.calculate_empty_slab',
            args=[self.calculator, input_data['slabParams']],
            inputs=['bulk_trajectory'],
            outputs=['slab_trajectory'],
        ))

        # Place Adsorbates
        tasks.append(fireworks.PyTask(
            func='catkit.flow.cohesive.calculate_adsorbates',
            args=[self.calculator,
                  input_data['slabParams'],
                  input_data['adsorbateParams'],
                  ],
            inputs=['bulk_trajectory'],
            outputs=['adsorbate_trajectories'],
        ))

        # Calculate Gas Phase Molecules
        tasks.append(fireworks.PyTask(
            func='catkit.flow.cohesive.calculate_molecules',
            args=[
                self.calculator,
                input_data['adsorbateParams'],
            ],
            inputs=[],
            outputs=['molecule_trajectories'],
        ))

        # Collected Energies and Report
        tasks.append(fireworks.PyTask(
            func='catkit.flow.cohesive.collect_energies',
            args=[input_data, ],
            inputs=[
                'bulk_trajectory',
                'slab_trajectory',
                'adsorbate_trajectories',
                'molecule_trajectories',
            ],
            outputs=[],
        ))

        if spec is None:
            spec = {'keys': keys, 'data': data}
        else:
            spec.update({'keys': keys, 'data': data})

        firework = fireworks.Firework(tasks,
                                      spec=spec,
                                      name=calculation_name,
                                      )

        workflow = fireworks.Workflow([firework])
        self.launchpad.add_wf(workflow)
