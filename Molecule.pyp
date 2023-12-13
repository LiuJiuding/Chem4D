################# Chem4D V1.0.0 ##################
import c4d
from PDT import PDTGeo, PDTFunction, C4DFunction, PDTTranslation
import os
import copy
from rdkit import Chem
import numpy as np

PLUGINID = 1059606
PLUGINNAME = "Chem4D Molecule"
VERSION = "1.0.0"
TITLE = f"{PLUGINNAME} v{VERSION}"
HELP = "molecular building tool"
# Dynamic group and parameters IDs
CHEM4D_DYNAMICGROUP = 1100
CHEM4D_DYNAMICGROUP_FIRSTPARAMETER = CHEM4D_DYNAMICGROUP + 1


class CHEM4DHelper(object):
    #read property file
    @staticmethod
    def ReadProperty(file):
        fprop = open(file,'r')
        dprop = {}
        for line in fprop:
            f2 = line.strip().split()
            dprop[f2[0]] = [int(f2[1]),np.array([int(f2[2]),int(f2[3]),int(f2[4])])/255]   #如果1个key有多个value，可以写成new_dict[f2[0]]  =  f2[1:]

        fprop.close()
        # print("read property ok")
        return dprop 
    
    @staticmethod
    def ReadMol(path):
        if path == '':
            return False
        if os.path.isfile(path) == False:
            print('path is not valid')
            return False
        #get molecular
        mol = Chem.MolFromMolFile(path,removeHs = False, strictParsing=False)
        if mol is None:
            mol = Chem.MolFromMol2File(path,removeHs = False)
            if mol is None:
                print('failed reading .mol/.mol2 file')
                return False
        else:
            print('successfully reading molecular file\n',path)
        return mol




class CHEM4D(c4d.plugins.ObjectData, CHEM4DHelper):

    def __init__(self):
        self.hasobj = False
        self.update = False
        self.read = False
        self.dprop = {}
        self.parameters = []
        self.symbols = []
        self.rootobj = c4d.BaseObject(5140)
        self.SetOptimizeCache(True)

    def Init(self, node):
        node[c4d.ID_MOL_ATOM_RADIUS_SCALE] = 1.0
        node[c4d.ID_MOL_ATOM_SEGMENTS] = 20
        node[c4d.ID_MOL_BOND_RADIUS] = 1.0
        node[c4d.ID_MOL_BUILD_MODE] = 0

        #read property file
        directory, _ = os.path.split(__file__)
        fn = os.path.join(directory, "elements.txt")
        self.dprop = self.ReadProperty(fn)
        return True

    def Read(self, node, hf, level):
        # Reads the number of dynamic parameters
        node[c4d.ID_MOL_ATOM_RADIUS_SCALE] = hf.ReadFloat32()
        node[c4d.ID_MOL_ATOM_SEGMENTS] = hf.ReadInt32()
        node[c4d.ID_MOL_BOND_RADIUS] = hf.ReadFloat32()
        self.hasobj = hf.ReadBool()
        self.read = True
        self.update = True
        count = hf.ReadInt32()
        # Reads the dynamic parameters value
        for idx in range(int(count/2)):
            value = hf.ReadFloat32()
            self.parameters.append(value)
            value = hf.ReadVector()
            self.parameters.append(value)
        return True

    def Write(self, node, hf):
        # Writes the number of dynamic parameters
        hf.WriteFloat32(node[c4d.ID_MOL_ATOM_RADIUS_SCALE])
        hf.WriteInt32(node[c4d.ID_MOL_ATOM_SEGMENTS])
        hf.WriteFloat32(node[c4d.ID_MOL_BOND_RADIUS])
        hf.WriteBool(self.hasobj)
        count = len(self.parameters)
        hf.WriteInt32(count)
        # Writes the dynamic parameters value
        for i in range(int(count/2)):
            hf.WriteFloat32(self.parameters[2*i])
            hf.WriteVector(self.parameters[2*i+1])
        return True

    def CopyTo(self, dest, snode, dnode, flags, trn):
        # Copies dynamic parameters value to the destination instance of ObjectData
        dest.update = copy.copy(self.update)
        dest.hasobj = copy.copy(self.hasobj)
        dest.read = copy.copy(self.read)
        dest.dprop = copy.copy(self.dprop)
        dest.symbols = copy.copy(self.symbols)
        dest.parameters = copy.copy(self.parameters)

        return True

    def GetDDescription(self, node, description, flags):
        # Loads the parameters from the description resource before adding dynamic parameters.
        if not description.LoadDescription(node.GetType()):
            return False

        # Get description single ID
        singleID = description.GetSingleDescID()

        # Declare dynamic group DescID
        dynamicGroupID = c4d.DescID(c4d.DescLevel(CHEM4D_DYNAMICGROUP, c4d.DTYPE_GROUP, node.GetType()))

        # Check if dynamic group needs to be added
        addDynamicGroup = singleID is None
        if not addDynamicGroup:
            addDynamicGroup = dynamicGroupID.IsPartOf(singleID)[0]

        # Adds dynamic group
        if addDynamicGroup:
            bc = c4d.GetCustomDataTypeDefault(c4d.DTYPE_GROUP)
            bc.SetString(c4d.DESC_NAME, "Atom Property")
            bc.SetInt32(c4d.DESC_COLUMNS, 1)
            if not description.SetParameter(dynamicGroupID, bc, c4d.DescID(c4d.DescLevel((c4d.ID_OBJECTPROPERTIES)))):
                return False

        # Declare REAL parameter container
        bc_real = c4d.GetCustomDataTypeDefault(c4d.DTYPE_REAL)
        bc_real.SetInt32(c4d.DESC_CUSTOMGUI, c4d.CUSTOMGUI_REAL)
        bc_real.SetFloat(c4d.DESC_MIN, 0.0)
        bc_real.SetFloat(c4d.DESC_MAX, 200.0)
        # bc_real.SetFloat(c4d.DESC_MINSLIDER, 0.0)
        # bc_real.SetFloat(c4d.DESC_MAXSLIDER, 200.0)
        bc_real.SetFloat(c4d.DESC_STEP, 1)
        bc_real.SetInt32(c4d.DESC_UNIT, c4d.DESC_UNIT_FLOAT)
        bc_real.SetInt32(c4d.DESC_ANIMATE, c4d.DESC_ANIMATE_ON)
        bc_real.SetBool(c4d.DESC_REMOVEABLE, False)

        bc_color = c4d.GetCustomDataTypeDefault(c4d.DTYPE_COLOR)
        bc_color.SetInt32(c4d.DESC_CUSTOMGUI, c4d.CUSTOMGUI_COLOR)

        # Initialize/Update parameters value list if needed
        parametersNum = len(self.symbols)*2
        parametersLen = len(self.parameters)

        if parametersNum == 0:
            self.parameters = []
        elif parametersLen != parametersNum:
            self.parameters.clear()
            for i in self.symbols:
                self.parameters.append(self.var_read_dprop[i][0])
                self.parameters.append(self.var_read_dprop[i][1])

        # Adds dynamic parameters
        idx = 0
        for symbol in self.symbols:
            descid = c4d.DescID(c4d.DescLevel(CHEM4D_DYNAMICGROUP_FIRSTPARAMETER+2*idx, c4d.DTYPE_REAL, node.GetType()))
            addParameter = singleID is None
            if not addParameter:
                addParameter = descid.IsPartOf(singleID)[0]

            if addParameter:
                name = "Atom Scale " + symbol
                bc_real.SetString(c4d.DESC_NAME, name)
                bc_real.SetString(c4d.DESC_SHORT_NAME, name)
                if not description.SetParameter(descid, bc_real, dynamicGroupID):
                    break

            descid = c4d.DescID(c4d.DescLevel(CHEM4D_DYNAMICGROUP_FIRSTPARAMETER+2*idx+1, c4d.DTYPE_COLOR, node.GetType()))
            addParameter = singleID is None
            if not addParameter:
                addParameter = descid.IsPartOf(singleID)[0]

            if addParameter:
                name = "Atom Color " + symbol
                bc_color.SetString(c4d.DESC_NAME, name)
                bc_color.SetString(c4d.DESC_SHORT_NAME, name)
                if not description.SetParameter(descid, bc_color, dynamicGroupID):
                    break
            idx = idx + 1

        # After dynamic parameters have been added successfully, return True and c4d.DESCFLAGS_DESC_LOADED with the input flags
        return True, flags | c4d.DESCFLAGS_DESC_LOADED

    def SetDParameter(self, node, id, data, flags):
        # Retrieves the parameter ID requested
        paramID = id[0].id

        # Retrieves the parameters count
        parametersLen = len(self.parameters)

        # Checks if passed parameter ID is a dynamic parameter
        if CHEM4D_DYNAMICGROUP_FIRSTPARAMETER <= paramID <= CHEM4D_DYNAMICGROUP+parametersLen:
            # Store the parameter data
            self.parameters[paramID-CHEM4D_DYNAMICGROUP_FIRSTPARAMETER] = data
            return True, flags | c4d.DESCFLAGS_SET_PARAM_SET

        return False

    def GetDParameter(self, node, id, flags):
        # Retrieves the parameter ID requested
        paramID = id[0].id

        # Retrieves the parameters count
        parametersLen = len(self.parameters)

        # Checks passed parameter ID is a dynamic parameter
        if CHEM4D_DYNAMICGROUP_FIRSTPARAMETER <= paramID <= CHEM4D_DYNAMICGROUP+parametersLen:
            # Retrieves the parameter data
            data = self.parameters[paramID-CHEM4D_DYNAMICGROUP_FIRSTPARAMETER]
            return True, data, flags | c4d.DESCFLAGS_GET_PARAM_GET
        return False

    def Message(self, node, type, data):
        if type==c4d.MSG_DESCRIPTION_COMMAND:
            if data['id'][0].id==c4d.ID_MOL_UPDATE:
                self.update = True
                node.SetDirty(c4d.DIRTYFLAGS_DATA)             
                return True
        return True

    def GetVirtualObjects(self, op, hh):
        #after click update button, cache all geometry under rootobj, change parameters of children only to speed up the plugin
        if self.hasobj == True and self.rootobj.GetDown() == None:
            self.update = True

        if self.update == True:
            scale = 100
            self.rootobj = c4d.BaseObject(5140)
            #check path and creat rdkit mol
            mol = self.ReadMol(op[c4d.ID_MOL_PATH])
            if mol == False:
                self.hasobj = False
                self.update = False
                return None


            #setup geo, geo_bond and self.parameters
            pos = mol.GetConformer().GetPositions() * scale
            geo = PDTGeo(pos)

            #atoms
            for i in range(geo.pointcount):
                symbol = mol.GetAtomWithIdx(i).GetSymbol()
                PDTFunction.setpointgroup(geo, symbol, i, 1)
                PDTFunction.setpointattrib(geo,'symbol', i, symbol)
                PDTFunction.setpointattrib(geo, "pscale", i, self.dprop[symbol][0])
                PDTFunction.setpointattrib(geo, "Cd", i, self.dprop[symbol][1])

            #bond
            geo_bond = copy.deepcopy(geo)
            for j in range(mol.GetNumBonds()):
                beginnum = mol.GetBondWithIdx(j).GetBeginAtom().GetIdx()
                endnum = mol.GetBondWithIdx(j).GetEndAtom().GetIdx()

                beginpos = geo_bond.GetP()[beginnum]
                endpos = geo_bond.GetP()[endnum]
                midpos = (beginpos + endpos) * 0.5

                midnum = PDTFunction.addpoint(geo_bond,midpos)

                bondbegin = PDTFunction.addprim(geo_bond,'polyline',[beginnum,midnum])
                bondend = PDTFunction.addprim(geo_bond,'polyline',[midnum,endnum])

                beginsymbol = mol.GetBondWithIdx(j).GetBeginAtom().GetSymbol()
                endsymbol = mol.GetBondWithIdx(j).GetEndAtom().GetSymbol()
                PDTFunction.setprimattrib(geo_bond,'symbol',bondbegin,beginsymbol)
                PDTFunction.setprimattrib(geo_bond,'Cd',bondbegin,self.dprop[beginsymbol][1])
                PDTFunction.setprimattrib(geo_bond,'symbol',bondend,endsymbol)
                PDTFunction.setprimattrib(geo_bond,'Cd',bondend,self.dprop[endsymbol][1])
            
            if self.read == False:
                self.parameters.clear()
                for i in geo.pointgroups:
                    self.parameters.append(self.dprop[i][0])
                    self.parameters.append(PDTTranslation.nptoc4d(self.dprop[i][1]))
                self.update = False
                self.hasobj = True
                
            self.symbols = geo.pointgroups
            #set null objects
            self.rootobj[c4d.ID_BASELIST_NAME] = os.path.basename(op[c4d.ID_MOL_PATH])

            rootatom = C4DFunction.CreateAtom(geo,"symbol")
            rootatom.InsertUnder(self.rootobj)

            sweep = C4DFunction.CreateBond(geo_bond,"symbol")
            sweep.SetName("bond")
            sweep.InsertUnder(self.rootobj)

            self.update = False
            self.hasobj = True
        
        if self.hasobj == False:
            return None
        #read static descriptions
        atomrad = op[c4d.ID_MOL_ATOM_RADIUS_SCALE]
        seg = op[c4d.ID_MOL_ATOM_SEGMENTS]
        bondrad = op[c4d.ID_MOL_BOND_RADIUS]

        #change parameters of cloner atom
        atomnull = self.rootobj.GetDownLast().GetDown().GetNext()
        clns = atomnull.GetChildren()
        # clns.reverse()
        for i in range(len(clns)):
            clns[i][c4d.ID_MG_TRANSFORM_COLOR] = self.parameters[2*i+1]
            clns[i].GetDown()[c4d.PRIM_SPHERE_SUB] = seg
            if op[c4d.ID_MOL_BUILD_MODE] == 2:
                clns[i].GetDown()[c4d.PRIM_SPHERE_RAD] = 10 * bondrad
            else:
                clns[i].GetDown()[c4d.PRIM_SPHERE_RAD] = self.parameters[2*i]*atomrad
        
        #change parameters of sweep bond
        bondnull = self.rootobj.GetDown()
        sweeps = bondnull.GetChildren()
        # sweeps.reverse()
        for i in range(len(sweeps)):
            sweeps[i].GetDown()[c4d.PRIM_CIRCLE_RADIUS] = 10 * bondrad
            sweeps[i][c4d.ID_BASEOBJECT_COLOR] = self.parameters[2*i+1]
        
        #always build bond, but show it if nessary
        if op[c4d.ID_MOL_BUILD_MODE] == 1:
            bondnull[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 1   #default is 2, hide is 1, show is 0
            bondnull[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 1
        else:
            bondnull[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 2
            bondnull[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 2
        return self.rootobj

# main
if __name__ == "__main__":
    # Retrieves the icon path
    directory, _ = os.path.split(__file__)
    fn = os.path.join(directory, "res", "omolreader.png")

    # Creates a BaseBitmap
    bmp = c4d.bitmaps.BaseBitmap()
    if bmp is None:
        raise MemoryError("Failed to create a BaseBitmap.")

    # Init the BaseBitmap with the icon
    if bmp.InitWith(fn)[0] != c4d.IMAGERESULT_OK:
        raise MemoryError("Failed to initialize the BaseBitmap.")

    # Registers the plugin
    c4d.plugins.RegisterObjectPlugin(id=PLUGINID,
                                      str=TITLE,
                                      g=CHEM4D,
                                      description="Omolreader", #file name without .res
                                      info=c4d.OBJECT_GENERATOR|c4d.OBJECT_USECACHECOLOR,
                                      icon=bmp)