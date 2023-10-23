from crystals import Crystal
path = r"G:\GitClone\Chem4D\cif\COD\LiCoO2.cif"
LCO = Crystal.from_cif(path)
print(LCO)
for atom in sorted(LCO.supercell(1,1,1)):
    print(repr(atom))
# print(LCO.__repr__())
