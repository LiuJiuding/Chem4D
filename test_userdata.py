import c4d
obj = doc.GetActiveObject()
bc_buildtype = c4d.GetCustomDataTypeDefault(c4d.DTYPE_LONG)

bc_buildtype[c4d.DESC_NAME] = "Build type"
bc_buildtype[c4d.DESC_CUSTOMGUI] = 200000281
bc = c4d.BaseContainer()
bc[0] = "a"
bc[1] = "b"
bc_buildtype[c4d.DESC_CYCLE] = bc
bc_buildtype[c4d.QUICKTAB_BAR] = True
bc_buildtype[c4d.QUICKTAB_BARTITLE] = "test"
# descId_buildtype = obj.AddUserData(bc_buildtype)
obj.SetUserDataContainer(c4d.DescID(c4d.DescLevel(c4d.ID_USERDATA),c4d.DescLevel(2)), bc_buildtype)
c4d.EventAdd()
