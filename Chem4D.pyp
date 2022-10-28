################# Chem4D V1.0.0 ##################
import c4d
from PDT import PDTGeo, PDTFunction, C4DFunction
import os
import copy
from rdkit import Chem

PLUGINID = 1059606
PLUGINNAME = "Chem4D"
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
            dprop[f2[0]] = [int(f2[1]),c4d.Vector(int(f2[2]),int(f2[3]),int(f2[4]))/255]   #如果1个key有多个value，可以写成new_dict[f2[0]]  =  f2[1:]

        fprop.close()

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
    @staticmethod
    def BuildPDTGeo(mol):
        pos = mol.GetConformer().GetPositions()
        largepos = []
        for p in pos:
            largepos.append(c4d.Vector(p[0]*100,p[1]*100,p[2]*100))
        pos = largepos
        geo = PDTGeo(pos)

        #atoms
        for i in range(geo.pointcount):
            PDTFunction.setpointgroup(geo, mol.GetAtomWithIdx(i).GetSymbol(), i, 1)
            PDTFunction.setpointattrib(geo,'symbol',i,mol.GetAtomWithIdx(i).GetSymbol())

        #bond
        geo_bond = geo
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
            PDTFunction.setprimattrib(geo_bond,'symbol',bondend,endsymbol)
        return geo,geo_bond

class CHEM4D(c4d.plugins.ObjectData, CHEM4DHelper):

    def __init__(self):
        self.hasobj = False
        self.update = False
        self.read = False
        self.read_dprop = {}
        self.parameters = []
        self.geo = PDTGeo()
        self.geo_bond = PDTGeo()
        self.rootobj = c4d.BaseObject(5140)
        self.SetOptimizeCache(True)

    def Init(self, node):
        node[c4d.ID_MOL_ATOM_RADIUS_SCALE] = 1.0
        node[c4d.ID_MOL_ATOM_SEGMENTS] = 20
        node[c4d.ID_MOL_BOND_RADIUS] = 1.0
        node[c4d.ID_MOL_BUILD_MODE] = 0

        #read property file
        directory, _ = os.path.split(__file__)
        fn = os.path.join(directory, "AtomProperties.txt")
        self.read_dprop = self.ReadProperty(fn)
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
        dest.read_dprop = copy.copy(self.read_dprop)
        dest.parameters = copy.copy(self.parameters)
        dest.geo = copy.copy(self.geo)
        dest.geo_bond = copy.copy(self.geo_bond)
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
        parametersNum = len(self.geo.pointgroups)*2
        parametersLen = len(self.parameters)

        if parametersNum == 0:
            self.parameters = []
        elif parametersLen != parametersNum:
            self.parameters.clear()
            for i in self.geo.pointgroups:
                self.parameters.append(self.var_read_dprop[i][0])
                self.parameters.append(self.var_read_dprop[i][1])

        # Adds dynamic parameters
        idx = 0
        for symbol in self.geo.pointgroups:
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
            if data['id'][0].id==c4d.ID_UPDATE:
                self.update = True
                node.SetDirty(c4d.DIRTYFLAGS_DATA)
                
                return True
        return True

    def GetVirtualObjects(self, op, hh):
        #after click update button, cache all geometry under rootobj, change parameters of children only to speed up the plugin
        if self.update == True:
            self.rootobj = c4d.BaseObject(5140)
            #check path and creat rdkit mol
            mol = self.ReadMol(op[c4d.ID_MOL_PATH])
            if mol == False:
                self.hasobj = False
                self.update = False
                return None

            #setup self.geo and self.parameters
            self.geo,self.geo_bond = self.BuildPDTGeo(mol)
            if self.read == False:
                self.parameters.clear()
                for i in self.geo.pointgroups:
                    self.parameters.append(self.read_dprop[i][0])
                    self.parameters.append(self.read_dprop[i][1])
                self.update = False
                self.hasobj = True
            
            #set null objects
            self.rootobj[c4d.ID_BASELIST_NAME] = os.path.basename(op[c4d.ID_MOL_PATH])
            rootatom = c4d.BaseObject(5140)
            rootatom[c4d.ID_BASELIST_NAME] = 'atom'
            rootbond = c4d.BaseObject(5140)
            rootbond[c4d.ID_BASELIST_NAME] = 'bond'
            rootatom.InsertUnder(self.rootobj)
            rootbond.InsertUnder(self.rootobj)

            #creat point cloud mesh
            pos = self.geo.GetP()
            pmesh = C4DFunction.postomesh(pos)
            pmesh[c4d.ID_BASELIST_NAME] = "PositionMesh"
            pmesh.InsertUnder(self.rootobj)

            #iterate atom symbols and build hierarchy
            indnum = 0
            rad = 1
            for i in self.geo.pointgroups:
                #creat atoms---------------------------------------------------------------------------------
                #creat point selection tags
                seltag = c4d.BaseTag(c4d.Tpointselection)
                seltag[c4d.ID_BASELIST_NAME] = i

                #the returned BaseSelect obj can be modified and automatically updated back to selection tag
                sel = seltag.GetBaseSelect()
                #select points from mark list, 1 is selected, 0 is unselected
                marklist = self.geo.GetPointGroup(i)
                sel.SetAll(marklist)
                
                #attach tag to mesh obj
                pmesh.InsertTag(seltag)

                #setup cloner
                sph = c4d.BaseObject(5160)
                sph[c4d.ID_BASELIST_NAME] = i
                sph[c4d.PRIM_SPHERE_RAD] = self.parameters[2*indnum]*rad
                sph[c4d.PRIM_SPHERE_TYPE] = 4
                sph[c4d.PRIM_SPHERE_SUB] = 10
                sphphone = c4d.BaseTag(5612)
                sphphone[c4d.PHONGTAG_PHONG_ANGLELIMIT] = 1
                sph.InsertTag(sphphone)
                cln = c4d.BaseObject(1018544)
                cln[c4d.ID_BASELIST_NAME] = i
                cln[c4d.MGCLONER_VOLUMEINSTANCES_MODE] = 1
                cln[c4d.ID_MG_TRANSFORM_COLOR] = self.parameters[2*indnum+1]
                cln[c4d.ID_MG_MOTIONGENERATOR_MODE] = 0
                cln[c4d.MG_OBJECT_LINK] = pmesh
                cln[c4d.MG_POLY_MODE_] = 0
                cln[c4d.MG_POLY_SELECTION] = i

                sph.InsertUnder(cln)
                cln.InsertUnder(rootatom)

                #creat bond splines------------------------------------------------------------------------
                #iterate bonds and get position list for atom type i
                #creat empty spline
                spl = c4d.BaseObject(5101)
                linepos = []
                for prim in self.geo_bond.primitives:
                    if prim.symbol == i:
                        for vertex in prim.vertices:
                            linepos.append(self.geo_bond.GetP()[vertex.pointnumber])
                #set point count and segment count
                spl.ResizeObject(len(linepos), int(len(linepos)/2))

                #set point positions and assign point to segment
                spl.SetAllPoints(linepos)
                for j in range(int(len(linepos)/2)):
                    spl.SetSegment(j, 2, False)
                
                #setup sweep
                spl[c4d.ID_BASELIST_NAME] = i
                spl.Message(c4d.MSG_UPDATE)
                circle = c4d.BaseObject(5181)
                circle[c4d.PRIM_CIRCLE_RADIUS] = 10
                sweep = c4d.BaseObject(5118)
                sweepphone = c4d.BaseTag(5612)
                sweepphone[c4d.PHONGTAG_PHONG_ANGLELIMIT] = 1
                sweep.InsertTag(sweepphone)
                sweep[c4d.ID_BASEOBJECT_USECOLOR] = 2
                sweep[c4d.ID_BASEOBJECT_COLOR] = self.parameters[2*indnum+1]
                sweep[c4d.ID_BASELIST_NAME] = i
                spl.InsertUnder(sweep)
                circle.InsertUnder(sweep)
                sweep.InsertUnder(rootbond)
                indnum = indnum +1
            self.update = False
            self.hasobj = True
            op.SetDirty(c4d.DIRTYFLAGS_DESCRIPTION)
        if self.hasobj == False:
            return None
        #read static descriptions
        atomrad = op[c4d.ID_MOL_ATOM_RADIUS_SCALE]
        seg = op[c4d.ID_MOL_ATOM_SEGMENTS]
        bondrad = op[c4d.ID_MOL_BOND_RADIUS]

        #change parameters of cloner atom
        clns = self.rootobj.GetDownLast().GetChildren()
        clns.reverse()
        for i in range(len(clns)):
            clns[i][c4d.ID_MG_TRANSFORM_COLOR] = self.parameters[2*i+1]
            clns[i].GetDown()[c4d.PRIM_SPHERE_SUB] = seg
            if op[c4d.ID_MOL_BUILD_MODE] == 2:
                clns[i].GetDown()[c4d.PRIM_SPHERE_RAD] = 10 * bondrad
            else:
                clns[i].GetDown()[c4d.PRIM_SPHERE_RAD] = self.parameters[2*i]*atomrad
        
        #change parameters of sweep bond
        sweeps = self.rootobj.GetDown().GetNext().GetChildren()
        sweeps.reverse()
        for i in range(len(sweeps)):
            sweeps[i].GetDown()[c4d.PRIM_CIRCLE_RADIUS] = 10 * bondrad
            sweeps[i][c4d.ID_BASEOBJECT_COLOR] = self.parameters[2*i+1]
        
        #always build bond, but show it if nessary
        if op[c4d.ID_MOL_BUILD_MODE] == 1:
            self.rootobj.GetDown().GetNext()[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 1   #default is 2, hide is 1, show is 0
            self.rootobj.GetDown().GetNext()[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 1
        else:
            self.rootobj.GetDown().GetNext()[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 2
            self.rootobj.GetDown().GetNext()[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 2
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