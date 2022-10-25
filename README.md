# Chem4D
 Use this python plugin to visualize moleculars in Cinema 4D, and creat interesting animation with MoGraph tools!
## Installation
1. Clone the repo and  put it in your C4D plugins folder "XXX\Maxon Cinema 4D R26\plugins"
2. Cut \_\_init__.py and PDT.py to C4D's python lib folder "C:\Users\XXX\AppData\Roaming\Maxon\Maxon Cinema 4D R26_042F1048\python39\libs"
3. Chem4D relies on third party python module to read mol files. Go to C4D installation folder "XXX\Maxon Cinema 4D R26\resource\modules\python\libs\python39.win64.framework", hold shift and click right mouse button to open powershell, type
`.\python -m pip install rdkit` and hit enter
4. Done!
## Change list
### v1.0.0
Features in initial release.<br>
1. Import .mol and .mol2 files
2. Change overall atom radius, atom segements and bond radius
3. Change atom color and radius for each element in the molecular
4. Chem4D reads the AtomProperties.txt file to set initial atom radius and color, you can add element symbols and modify this file. e.g. H radius R G B
