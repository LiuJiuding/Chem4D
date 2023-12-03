import c4d
obj = c4d.BaseObject(1061833)
obj[c4d.ID_CIF_PATH] = r"G:\GitClone\Chem4D\sample\LiFePO4.cif"
doc.InsertObject(obj)
c4d.EventAdd()