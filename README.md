# Chem4D
 Chem4D is a python plugin collection for Cinema 4D, which enables you to import chemical structure files in Cinema 4D.
## Installation
1. Download the release package and unzip it in your C4D plugins folder "XXX\Maxon Cinema 4D RXX\plugins"
2. Cut \_\_init__.py and PDT.py to C4D's python lib folder "C:\Users\XXX\AppData\Roaming\Maxon\Maxon Cinema 4D RXX_XXXXXX\python39\libs"
3. Chem4D relies on third party python module `rdkit` and `pymatgen`. To install these to C4D python environment, go to C4D installation folder "XXX\Maxon Cinema 4D R26\resource\modules\python\libs\python39.win64.framework" or "XXX\Maxon Cinema 4D 2024\resource\modules\python\libs\win64", backup this folder, then hold shift and click right mouse button to open powershell, input following command and hit enter one by one.<br>
`./python -m ensurepip --upgrade`<br>
`./python -m pip install rdkit`<br>
`./python -m pip install pymatgen`
4. Done!
## Tips
* Use .mol2 file for `Molecule` whenever possible, it supports more molecular data.
* Press C (make editable) and you will find that this molecule is built with Cloners.
* Open C4D console to see outputs of this plugin when it doesn't work right.
* It is a good habit to backup your work when using Chem4D, it may be unstable :)
* The plugin is only tested on C4D R26.107, other versions may work.
## Change list
### v1.0.0
Features in initial release.<br>
1. Import .mol and .mol2 files.
2. Change overall atom radius, atom segements and bond radius.
3. Change atom color and radius for each element in the molecular.
4. Chem4D reads the AtomProperties.txt file to set initial atom radius and color, you can add element symbols and modify this file. e.g. H radius R G B

Issues.
* Can NOT undo if you delete the Chem4D object in Attribute Manager.
### v2.0.0
Major feature release.<br>
Add a new plugin `Crystal` for reading small structure cif file. Original mol reader plugin has been renamed as `Molecule`.
#### Features
1. `Crystal` Import inorganic structure .cif file.
2. `Crystal` Change overall atom radius, atom segements and bond radius.
3. `Crystal` Change atom color and radius for each element in the crystal.
4. `Crystal` and `Molecule` Reads the elements.txt file to set initial atom radius and color, you can add element symbols and modify this file. e.g. H radius R G B
5. `Crystal` You can choose where to generate polyhedral. <br>
Guess mode will guess charge on each site, and polyhedral will be build on site with positive charge. Do not use guess mode when a unit cell is big, it is very slow e.g. ZIF-8<br>
Symbol or Label mode will use input symbol or label as polyhedral center. Label can be found in cif file under `_atom_site_label`, and symbol is pure element symbol. Multiple words should be seperated by space.<br>
6. `Crystal` Change cell a, b, c dimention
7. `Molecule` Adopted to new PDT.py
#### Limitations
1. `Crystal` Cannot identify molecule piece on the cell boundary. To get full molecule, generate `3*3*3` structure and delete outer atoms.
2. `Crystal` Due to the brute force algorithm to remove overlaped atoms and bonds, it is slow to generate supercell larger than 100.
3. `Crystal` `_atom_site_occupancy` must be 1.
4. `Crystal` Some unit cell structure may missing a corner, this is a bug in pymatgen. e.g. Fe3O4_spinel.cif
5. `Crystal` and `Molecule` Do not use Enabel toggle in ObjectManager, there is a bug that has not been resolved.
6. `Crystal` and `Molecule` Dynamic parameters (Atom property) will loss when undo during regenrating geometry


