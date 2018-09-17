## Using the cathub cli

Run `cathub`, like so

    cathub --help

or with any of its sub-commands, like so

    cathub make_folders_template --help

## Examples

To organize DFT Input files into a catalysis-hub.org compliant folder structure

    cathub organize folder_name -d gas_molecules_folder -a CO2,CO,C -c VASP-5.4.4 -x PBE+U -f 111 -S Ni_fcc -v -k -r C,O2 

To create an .json input file

    cathub make_folders_template project1.json --create-template

To create a folder structures from a .json input file

    cathub make_folders_template project1.json

Querying the Catalysis Hub database:

    cathub reactions -q reactants=CO -q chemicalComposition=~Pt

    cathub publications -q title=~Evolution -q year=2017

Reading folders into sqlite3 db file:

    cathub folder2db <foldername>

Sending the data to the Catalysis Hub server:

    cathub db2server <dbfile>
