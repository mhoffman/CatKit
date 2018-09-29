from netrc import netrc
import copy
import importlib
import json
import pprint
import subprocess
import tempfile
import os
import zipfile
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
import catkit.hub.organize

USERHANDLE = 'mjhoffmann@gmail.com'
DFT_CODE = 'QE-5.1'
## DEFAULT PARAMETERS FOR CALCULATOR
## CAN BE OVERRIDEN BY CALCULATOR_PARAMETERS
## ON INPUT DATA
BULK_CALCULATOR = {
    'input_dft': 'BEEF-vdw',    # Requires Johnannes' version of QE
    'calculation': 'vc-relax',
    'kspacing': 0.02,
    'forc_conv_thr': 0.05,      # Force convergence (same units as ASE)
    'nbnd': -10,                # A safety precausion
    'ecutwfc': 800,             # Very high, recommended by Martin
    'ecutrho': 8000,
    'conv_thr': 1e-4,           # More reasonable considering expected DFT accuracy
    'nosym': True
    }

SLAB_CALCULATOR = {
    'input_dft': 'BEEF-vdw',
    'calculation': 'relax',
    'kspacing': 0.04,
    'forc_conv_thr': 0.05,
    'nbnd': -10,
    'ecutwfc': 500,
    'ecutrho': 5000,
    'conv_thr': 1e-4,
    'degauss': 0.2,
    'mixing_beta': 0.2,
    'dipfield': True,
    'nosym': True
        }

MOLECULE_CALCULATOR = {
    'input_dft': 'BEEF-vdw',
    'calculation': 'relax',
    'forc_conv_thr': 0.05,
    'ecutwfc': 500,
    'ecutrho': 5000,
    'conv_thr': 1e-4,
    'degauss': 0.01,
    'nosym': True
        }


def open_string(string, format='cif'):
    with StringIO.StringIO() as temp_file:
        temp_file.write(string)
        temp_file.seek(0)
        atoms = ase.io.read(
            temp_file,
            format=format,
        )
    return atoms

def to_string(atoms, format='json'):
    with StringIO.BytesIO() as temp_file:
        ase.io.write(temp_file, atoms, format=format)
        temp_file.seek(0)
        res = temp_file.getvalue()
    return res



def str_to_class(full_classname):
    mod_name, classname = full_classname.rsplit('.', 1)
    return getattr(importlib.import_module(mod_name),
                   classname)


def setup_calculator(
    atoms,
    calculator_name
):
    # Setting up the calculator
    calculator = str_to_class(calculator_name)
    if calculator_name == 'espresso.Espresso':
        calc = calculator(atoms, **atoms.info.get('calculator_parameters', {}))

    else:
        calc = calculator(**atoms.info.get('calculator_parameters', {}))
        atoms.set_calculator(calc)


def calculate_bulk(
        calculator,
        dft_params,
        bulk_params,
):
    print("============================")
    print("IN CALCULATE BULK")
    print(calculator)
    print(bulk_params)

    atoms = open_string(bulk_params['wyckoff']['cif'])
    atoms.info.setdefault('calculator_parameters', {}) \
            .update(BULK_CALCULATOR)
    atoms.info.setdefault('calculator_parameters', {}) \
            .update(dft_params)

    setup_calculator(atoms, calculator)
    energy = atoms.get_potential_energy()

    print(energy)
    print("============================")
    atoms.info.update({
        
        })
    return fwio.atoms_to_encode(atoms)


def calculate_empty_slab(
        calculator,
        dft_params,
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
    slab.info.setdefault('calculator_parameters', {}) \
            .update(SLAB_CALCULATOR)
    slab.info.setdefault('calculator_parameters', {}) \
            .update(dft_params)

    slab.info.update(slab_params)

    setup_calculator(slab, calculator)
    slab.get_potential_energy()

    print("============================")
    return fwio.atoms_to_encode(slab)


def calculate_adsorbates(
    calculator,
    dft_params,
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

    print("CREATED " + str(len(slabs)) + " SLABS")
    pprint.pprint(slabs)
    print("FOUND " + str(len(adsorbate_params['equations'])) + " EQUATIONS")
    for i, slab in enumerate(slabs):
        setup_calculator(slab, calculator)
        slab.get_potential_energy()
        slab.info.update({
            'equation': adsorbate_params['equations'][i],
            'adsorbate': adsorbate_params['adsorbate'],
            })

    encoded_slabs = [
        fwio.atoms_to_encode(slab)
        for slab in slabs
    ]
    for i, slab in enumerate(encoded_slabs):
        print(" ===> SLAB")
        print(adsorbate_params['equations'][i])
        print(slab)
    print("============================")

    return json.dumps([
        fwio.atoms_to_encode(slab)
        for slab in slabs
    ])


def calculate_molecules(
        calculator,
        dft_params,
        adsorbate_params,
):
    print("============================")
    print("IN CALCULATE MOLECULES")
    # print(calculator)
    #print(adsorbate_params)
    input_format = adsorbate_params['format']
    molecules = copy.deepcopy(adsorbate_params['molecules'])

    trajectories = {}

    for name, input_string in molecules.items():
        molecule = open_string(input_string, input_format)
        setup_calculator(molecule, calculator)
        molecule.get_potential_energy()
        trajectories[name] = fwio.atoms_to_encode(molecule)


    print("============================")

    return json.dumps(trajectories)


def collect_energies(
    input_data,
    dft_params,
    bulk_trajectory,
    slab_trajectory,
    adsorbate_trajectories,
    molecule_trajectories,
):
    # Deserialize JSON strings
    bulk_trajectory = json.loads(bulk_trajectory)
    slab_trajectory = json.loads(slab_trajectory)
    adsorbate_trajectories =  [
            json.loads(x)
            for x in json.loads(adsorbate_trajectories)
            ]
    molecule_trajectories = json.loads(molecule_trajectories)


    print("============================")
    print("IN COLLECTING ENERGIES")
    #print(input_data)
    #print("## BULK")
    #print(bulk_trajectory)
    #print(type(bulk_trajectory))
    print("## EMPTY SLAB")
    print(slab_trajectory)
    print(type(slab_trajectory))
    #print("## ADSORBATES")
    #print(adsorbate_trajectories)
    #print(type(adsorbate_trajectories))
    #print("## MOLECULES")
    #print(molecule_trajectories)


    with open('collected_energies.json', 'w') as outfile:
        outfile.write(json.dumps({
            'input_data': input_data,
            'bulk_trajectory': bulk_trajectory,
            'slab_trajectory': slab_trajectory,
            'adsorbate_trajectories': adsorbate_trajectories,
            'molecule_trajectories': molecule_trajectories,
        }))
        # outfile.write(json.dumps(kwargs))

    ##########################################################
    # Rationale:
    # In order to use the cathub ... toolchain
    # we first create a zip file in memory
    # since it offers a simple API for creating deeply nested
    # directory structures
    # then we unzip to a temporary directory
    # and use the cathub CLI to upload the result
    ##########################################################
    zip_mem_file = StringIO.BytesIO()
    zf = zipfile.ZipFile(zip_mem_file, 'w')

    publication_path = ('calculation/'  \
                + '/publication.txt'
                ).format(**locals())

    zf.writestr(publication_path, catkit.hub.organize.PUBLICATION_TEMPLATE)



    print(bulk_trajectory)
    dft_functional = bulk_trajectory['info']['calculator_parameters']['input_dft']
    dft_code = DFT_CODE

    composition = ''.join(input_data['bulkParams']['elements'])
    structure = input_data['bulkParams']['wyckoff']['name']

    # save bulk structures
    bulk_path = ('calculation/{dft_code}'  \
                + '/{dft_functional}' \
                + '/{composition}__{structure}'
                + '/bulk' \
                + '/bulk.json'
                ).format(**locals())
    zf.writestr(bulk_path, json.dumps(bulk_trajectory))

    # save gas phase molecules
    for molecule, molecule_structure in molecule_trajectories.items():
        molecule_path = ('calculation/{dft_code}'  \
                    + '/{dft_functional}' \
                    + '/gas' \
                    + '/{molecule}.json'
                    ).format(**locals())
        print("RAW MOLECULE STRUCTURE")
        print(molecule_structure)
        print(type(molecule_structure))
        print('#######')
        molecule = json.loads(molecule_structure)['trajectory']
        molecule = fwio.encode_to_atoms(molecule_structure)
        print("MOLECULE")
        print(molecule)
        print(type(molecule))
        print("=====================")
        zf.writestr(
                molecule_path,
                to_string(molecule),
                )

    millerX = input_data['slabParams']['millerX']
    millerY = input_data['slabParams']['millerY']
    millerZ = input_data['slabParams']['millerZ']

    facet = '{millerX}{millerY}{millerZ}'.format(**locals())

    # save empty slab
    slab_path = ('calculation/{dft_code}'  \
                + '/{dft_functional}' \
                + '/{composition}__{structure}'
                + '/{facet}' \
                + '/empty_slab.json' 
                ).format(**locals())
    slab = fwio.encode_to_atoms(json.dumps(slab_trajectory))
    zf.writestr(slab_path, to_string(slab))

    # save adsorbate structure
    for adsorbate in adsorbate_trajectories:
        equation = adsorbate['info']['equation']
        adsorbate_name = adsorbate['info']['adsorbate']
        adsorbate_path = ('calculation/{dft_code}'  \
                    + '/{dft_functional}' \
                    + '/{composition}__{structure}'
                    + '/{facet}' \
                    + '/{equation}' \
                    + '/{adsorbate_name}.json' 
                    ).format(**locals())
        print("ADSORBATE")
        print(adsorbate)
        print(type(adsorbate))
        adsorbate_atoms = fwio.encode_to_atoms(json.dumps(adsorbate))
        print("===========================")
        zf.writestr(adsorbate_path, to_string(adsorbate_atoms))






    temp_dir = tempfile.mkdtemp()
    zf.extractall(temp_dir)
    zf.close()

    pl = subprocess.Popen((
        'tree -hs {temp_dir}'.format(**locals())
        ).split(),
        stdout=subprocess.PIPE,
        )
    print(pl.communicate()[0])


    userhandle = USERHANDLE
    print(
        'cathub folder2db --energy-limit 1000 --userhandle {userhandle} {temp_dir}'.format(**locals())
            )
    p1 = subprocess.Popen((
        'cathub folder2db --energy-limit 1000 --userhandle {userhandle} {temp_dir}'.format(**locals())
            ).split(),
        stdout=subprocess.PIPE,
            )
    print(p1.communicate()[0])

    db_file = os.path.join(*[
        temp_dir,
        'DoeTest2018.db'

        ])
    pl = subprocess.Popen((
        'tree -hs {temp_dir}'.format(**locals())
        ).split(),
        stdout=subprocess.PIPE,
        )
    print(pl.communicate()[0])



    print(
        'cathub db2server {db_file}'.format(**locals())
            )
    p2 = subprocess.Popen((
        'cathub db2server {db_file}'.format(**locals())
        ).split(),
        stdout=subprocess.PIPE,
        )
    print(p2.communicate()[0])

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
                input_data['dftParams'],
                input_data['bulkParams'],
            ],
            outputs=['bulk_trajectory'],
            auto_kwargs=True,
        ))

        # Cut Empty Slab
        tasks.append(fireworks.PyTask(
            func='catkit.flow.cohesive.calculate_empty_slab',
            args=[
                self.calculator,
                input_data['dftParams'],
                input_data['slabParams'],
                ],
            inputs=['bulk_trajectory'],
            outputs=['slab_trajectory'],
        ))

        # Place Adsorbates
        tasks.append(fireworks.PyTask(
            func='catkit.flow.cohesive.calculate_adsorbates',
            args=[self.calculator,
                  input_data['dftParams'],
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
                input_data['dftParams'],
                input_data['adsorbateParams'],
            ],
            inputs=[],
            outputs=['molecule_trajectories'],
        ))

        # Collected Energies and Report
        tasks.append(fireworks.PyTask(
            func='catkit.flow.cohesive.collect_energies',
            args=[
                input_data,
                input_data['dftParams'],
                ],
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
