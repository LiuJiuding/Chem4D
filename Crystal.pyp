################# Chem4D Crystal V1.0.0 ##################
from typing import Optional
import c4d
from c4d import utils
import os
import copy

import PDT
from PDT import PDTGeo, PDTFunction, C4DFunction, PDTTranslation

import scipy.spatial as spt
import pymatgen.core as pmg
import pymatgen.analysis.local_env as pmg_env
import numpy as np

doc: c4d.documents.BaseDocument  # The active document
op: Optional[c4d.BaseObject]  # The active object, None if unselected

PLUGINID = 	1061833
PLUGINNAME = "Chem4D Crystal"
VERSION = "1.0.0"
TITLE = f"{PLUGINNAME} v{VERSION}"
HELP = "crystal building tool"
# Dynamic group and parameters IDs
CHEM4D_DYNAMICGROUP = 1100
CHEM4D_DYNAMICGROUP_FIRSTPARAMETER = CHEM4D_DYNAMICGROUP + 1
        
class CHEM4DHelper(object):

    @staticmethod
    def GetSymbol(site):
        symbol = ""
        if type(site.specie)==pmg.periodic_table.Species:
            symbol = str(site.specie.element)
        else:
            symbol = str(site.specie)
        return symbol
    
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
    def ConvexHull(convex_pos:list[c4d.Vector]) -> c4d.PolygonObject:
        makeconvex = True
        try:
            hull = spt.ConvexHull(PDTTranslation.c4dtonp(convex_pos))
        except spt.QhullError:
            makeconvex = False
        if makeconvex:
        # pos = nptoc4d(convex_pos)
            poly = c4d.PolygonObject(0,0)
            poly.ResizeObject(len(convex_pos),len(hull.simplices))
            poly.SetAllPoints(convex_pos)
            i = 0
            for simplex in hull.simplices:
                ptnum = [simplex[0],simplex[1],simplex[2]]
                a = convex_pos[simplex[1]] - convex_pos[simplex[0]]
                b = convex_pos[simplex[2]] - convex_pos[simplex[0]]
                normal = a.Cross(b)
                angle = 0
                for j in range(len(convex_pos)):
                    if j not in ptnum:
                        v = convex_pos[j] - convex_pos[simplex[0]]
                        angle = utils.RadToDeg(utils.VectorAngle(v,normal))
                        break
                if angle > 90:
                    poly.SetPolygon(i,c4d.CPolygon(simplex[0],simplex[1],simplex[2]))
                else: poly.SetPolygon(i,c4d.CPolygon(simplex[0],simplex[2],simplex[1]))
                i+=1
            poly.Message(c4d.MSG_UPDATE)
            return poly
        else:
            print("skip building convexhull ")
            return c4d.PolygonObject(0,0)
    
    @staticmethod
    def ReadCif(path):
            if path == '':
                return False
            if os.path.isfile(path) == False:
                print('path is not valid')
                return False
            cryst = pmg.Structure.from_file(path)
            if cryst is None or False:
                print('failed reading .cif file')
                return False
            else:
                print('successfully reading cif file\n',path)
            return cryst

class CHEM4DCrystal(c4d.plugins.ObjectData, CHEM4DHelper):

    def __init__(self):
        self.hasobj = False
        self.update = False
        self.read = False
        self.dprop = {}
        self.parameters = []
        self.symbols = []
        self.label_dict = {}
        self.rootobj = c4d.BaseObject(5140)
        self.SetOptimizeCache(True)

    def Init(self, node):
        node[c4d.ID_CIF_ATOM_RADIUS_SCALE] = 1.0
        node[c4d.ID_CIF_ATOM_SEGMENTS] = 20
        node[c4d.ID_CIF_BOND_RADIUS] = 1.0
        node[c4d.ID_CIF_BUILD_MODE] = 0
        node[c4d.ID_CIF_SHOW_CELL] = 1
        node[c4d.ID_CIF_SHOW_BOND] = 1
        node[c4d.ID_CIF_SHOW_POLY] = 1
        node[c4d.ID_CIF_A] = 1
        node[c4d.ID_CIF_B] = 1
        node[c4d.ID_CIF_C] = 1
        node[c4d.ID_CIF_CENTER_TYPE] = 0

        #read property file
        directory, _ = os.path.split(__file__)
        fn = os.path.join(directory, "elements.txt")
        self.dprop = self.ReadProperty(fn)
        return True

    def Read(self, node, hf, level):
        # Reads the number of dynamic parameters
        node[c4d.ID_CIF_ATOM_RADIUS_SCALE] = hf.ReadFloat32()
        node[c4d.ID_CIF_ATOM_SEGMENTS] = hf.ReadInt32()
        node[c4d.ID_CIF_BOND_RADIUS] = hf.ReadFloat32()
        node[c4d.ID_CIF_SHOW_CELL] = hf.ReadBool()
        node[c4d.ID_CIF_SHOW_BOND] = hf.ReadBool()
        node[c4d.ID_CIF_SHOW_POLY] = hf.ReadBool()
        node[c4d.ID_CIF_BUILD_MODE] = hf.ReadInt32()
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
        hf.WriteFloat32(node[c4d.ID_CIF_ATOM_RADIUS_SCALE])
        hf.WriteInt32(node[c4d.ID_CIF_ATOM_SEGMENTS])
        hf.WriteFloat32(node[c4d.ID_CIF_BOND_RADIUS])
        hf.WriteBool(node[c4d.ID_CIF_SHOW_CELL])
        hf.WriteBool(node[c4d.ID_CIF_SHOW_BOND])
        hf.WriteBool(node[c4d.ID_CIF_SHOW_POLY])
        hf.WriteInt32(node[c4d.ID_CIF_BUILD_MODE])
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
        dest.label_dict = copy.copy(self.label_dict)
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
                self.parameters.append(self.dprop[i][0])
                self.parameters.append(PDTTranslation.nptoc4d(self.dprop[i][1]))

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
            if data['id'][0].id==c4d.ID_CIF_UPDATE:
                # print("click updata")
                self.update = True
                node.SetDirty(c4d.DIRTYFLAGS_DATA)             
                return True
        return True

    def GetVirtualObjects(self, op, hh):
        # after click update button, cache all geometry under rootobj, change parameters of children only to speed up the plugin
        if self.hasobj == True and self.rootobj.GetDown() == None:
            self.update = True

        if self.update == True:
            self.rootobj = c4d.BaseObject(5140)
            scale = np.array([100.0,100.0,-100.0])
            delta = 0.0001

            # read cif, get dictionary and label list----------------------------------------
            cryst = self.ReadCif(op[c4d.ID_CIF_PATH])
            if cryst == False:
                return None
            
            charge_dict = {}

            
            sym = []
            if op[c4d.ID_CIF_CENTER_TYPE] == 0:
                cryst.add_oxidation_state_by_guess()
                for site in cryst:
                    charge_dict[site.label] = str(site.specie)[-1]
                    self.label_dict[site.label] = str(site.specie.element)
                    sym.append(str(site.specie.element))
            elif op[c4d.ID_CIF_CENTER_TYPE] == 1:
                poly_center = op[c4d.ID_CIF_POLY_CENTER].split()
                for site in cryst:
                    if self.GetSymbol(site) in poly_center:
                        charge_dict[site.label] = "+"
                    else: charge_dict[site.label] = "-"
                    self.label_dict[site.label] = self.GetSymbol(site)
                    sym.append(self.GetSymbol(site))
            elif op[c4d.ID_CIF_CENTER_TYPE] == 2:
                poly_center = op[c4d.ID_CIF_POLY_CENTER].split()
                for site in cryst:
                    if site.label in poly_center:
                        charge_dict[site.label] = "+"
                    else: charge_dict[site.label] = "-"
                    self.label_dict[site.label] = self.GetSymbol(site)
                    sym.append(self.GetSymbol(site))
            self.symbols = list(set(sym))

            if self.read == False:
                self.parameters.clear()
                for i in self.symbols:
                    self.parameters.append(self.dprop[i][0])
                    self.parameters.append(PDTTranslation.nptoc4d(self.dprop[i][1]))
                self.update = False
                self.hasobj = True

            # print(charge_dict,self.label_dict)
            anion_label = []
            cation_label = []
            for key,value in charge_dict.items():
                if value == "-":
                    anion_label.append(key)
                if value == "+":
                    cation_label.append(key)

            # set null objects----------------------------------------------------------------

            self.rootobj[c4d.ID_BASELIST_NAME] = os.path.basename(op[c4d.ID_CIF_PATH])
            rootpoly = c4d.BaseObject(5140)
            rootpoly[c4d.ID_BASELIST_NAME] = 'poly'
            rootpoly.InsertUnder(self.rootobj)

            # setup unit cell line-------------------------------------------------------------------
            geo_cell = PDTGeo()
            axises = ["a","b","c"]
            colors = [np.array([1.0,0.0,0.0]),np.array([0.0,1.0,0.0]),np.array([0.0,0.0,1.0])]
            for i in range(3):
                beginpos = np.array([0.0,0.0,0.0])
                endpos = cryst.lattice.matrix[i]*scale
                beginnum = PDTFunction.addpoint(geo_cell,beginpos)
                endnum = PDTFunction.addpoint(geo_cell,endpos)
                celledge = PDTFunction.addprim(geo_cell,'polyline',[beginnum,endnum])
                PDTFunction.setprimattrib(geo_cell,'label',celledge,axises[i])
                PDTFunction.setprimattrib(geo_cell,'Cd',celledge,colors[i])

                beginpos = np.array([0.0,0.0,0.0]) + cryst.lattice.matrix[(i+1)%3]*scale
                endpos = cryst.lattice.matrix[i]*scale + cryst.lattice.matrix[(i+1)%3]*scale
                beginnum = PDTFunction.addpoint(geo_cell,beginpos)
                endnum = PDTFunction.addpoint(geo_cell,endpos)
                celledge = PDTFunction.addprim(geo_cell,'polyline',[beginnum,endnum])
                PDTFunction.setprimattrib(geo_cell,'label',celledge,axises[i])
                PDTFunction.setprimattrib(geo_cell,'Cd',celledge,colors[i])

                beginpos = np.array([0.0,0.0,0.0]) + cryst.lattice.matrix[(i+2)%3]*scale
                endpos = cryst.lattice.matrix[i]*scale + cryst.lattice.matrix[(i+2)%3]*scale
                beginnum = PDTFunction.addpoint(geo_cell,beginpos)
                endnum = PDTFunction.addpoint(geo_cell,endpos)
                celledge = PDTFunction.addprim(geo_cell,'polyline',[beginnum,endnum])
                PDTFunction.setprimattrib(geo_cell,'label',celledge,axises[i])
                PDTFunction.setprimattrib(geo_cell,'Cd',celledge,colors[i])

                beginpos = np.array([0.0,0.0,0.0]) + cryst.lattice.matrix[(i+1)%3]*scale + cryst.lattice.matrix[(i+2)%3]*scale
                endpos = cryst.lattice.matrix[i]*scale + cryst.lattice.matrix[(i+1)%3]*scale + cryst.lattice.matrix[(i+2)%3]*scale
                beginnum = PDTFunction.addpoint(geo_cell,beginpos)
                endnum = PDTFunction.addpoint(geo_cell,endpos)
                celledge = PDTFunction.addprim(geo_cell,'polyline',[beginnum,endnum])
                PDTFunction.setprimattrib(geo_cell,'label',celledge,axises[i])
                PDTFunction.setprimattrib(geo_cell,'Cd',celledge,colors[i])
            
            # make super cell to fill unit cell----------------------------------------------------
            cryst.make_supercell([2,2,2])
            m = cryst.lattice.matrix
            geo = PDTGeo()
            print("initialize NN")
            NN = pmg_env.CrystalNN(distance_cutoffs=None,x_diff_weight=0,porous_adjustment=True)

            # build atom inside boundary-----------------------------------------------------------
            print("find inside")
            inside = []
            for i in range(len(cryst)):
                f = cryst[i].frac_coords
                if f[0] <= 0.5+delta and f[1] <= 0.5+delta and f[2] <= 0.5+delta:
                    inside.append(i)
                    p = cryst[i].coords * scale
                    ptnum = PDTFunction.addpoint(geo, p)
                    symbol = self.GetSymbol(cryst[i])
                    PDTFunction.setpointattrib(geo, "symbol", ptnum, symbol)
                    PDTFunction.setpointattrib(geo, "label", ptnum, cryst[i].label)
                    PDTFunction.setpointattrib(geo, "Cd", ptnum, self.dprop[symbol][1])
                    PDTFunction.setpointattrib(geo, "pscale", ptnum, self.dprop[symbol][0])
                    PDTFunction.setpointgroup(geo, cryst[i].label, ptnum, 1)
            labels = geo.pointgroups

            # build poly and bond-------------------------------------------------------------------
            geo_bond = PDTGeo()

            poly_dict = {}
            for label in cation_label:
                rootconvex = c4d.BaseObject(5140)
                rootconvex[c4d.ID_BASELIST_NAME] = label
                rootconvex.InsertUnder(rootpoly)
                poly_dict[label] = rootconvex

            bondarr = np.empty(shape=(0,8))
            for i in inside:
                print("iterate site: "+str(i))
                dicts = NN.get_nn_info(cryst,i)
                # print("get NN info: "+str(i))
                # begin property
                beginpos = cryst[i].coords * scale
                beginlabel = cryst[i].label
                beginsymbol = self.GetSymbol(cryst[i])
                pos = []
                for dic in dicts:

                    # end property
                    endpos = dic["site"].coords * scale
                    pos.append(endpos)
                    endlabel = dic["site"].label
                    endsymbol = self.GetSymbol(dic["site"])

                    midpos = (beginpos + endpos) * 0.5

                    p = dic["site"].coords
                    f = np.linalg.solve(m.T,p)
                    
                    outside = f[0] > 0.5+delta or f[1] > 0.5+delta or f[2] > 0.5+delta or f[0] < 0.0 or f[1] < 0.0 or f[2] < 0.0
                    # print(f,outside)
                    # 中心是阴离子且配体在外，不创建配体原子和键
                    if not(charge_dict[cryst[i].label] == "-" and outside):
                        # 中心是阳离子且配体在外，创建配体原子
                        if charge_dict[cryst[i].label] == "+" and outside:
                            # atom
                            ptnum = PDTFunction.addpoint(geo, endpos)
                            PDTFunction.setpointattrib(geo, "symbol", ptnum, endsymbol)
                            PDTFunction.setpointattrib(geo, "label", ptnum, endlabel)
                            PDTFunction.setpointattrib(geo, "Cd", ptnum, self.dprop[endsymbol][1])
                            PDTFunction.setpointattrib(geo, "pscale", ptnum, self.dprop[symbol][0])
                            PDTFunction.setpointgroup(geo, endlabel, ptnum, 1)
                        # bond
                        beginnum = PDTFunction.addpoint(geo_bond,beginpos)
                        midnum = PDTFunction.addpoint(geo_bond,midpos)
                        endnum = PDTFunction.addpoint(geo_bond,endpos)

                        bondbegin = PDTFunction.addprim(geo_bond,'polyline',[beginnum,midnum])
                        PDTFunction.setprimattrib(geo_bond,'label',bondbegin,beginlabel)
                        PDTFunction.setprimattrib(geo_bond,'Cd',bondbegin,self.dprop[beginsymbol][1])

                        bondend = PDTFunction.addprim(geo_bond,'polyline',[midnum,endnum])
                        PDTFunction.setprimattrib(geo_bond,'label',bondend,endlabel)
                        PDTFunction.setprimattrib(geo_bond,'Cd',bondend,self.dprop[endsymbol][1])

                # poly
                pos = PDTTranslation.listtoc4d(pos)
                if charge_dict[cryst[i].label] == "+" and len(pos)>3:
                    poly = self.ConvexHull(pos)
                    poly.SetName(cryst[i].label)
                    poly[c4d.ID_BASEOBJECT_USECOLOR] = 2
                    poly[c4d.ID_BASEOBJECT_COLOR] = PDTTranslation.nptoc4d(self.dprop[self.GetSymbol(cryst[i])][1])
                    poly[c4d.ID_BASEOBJECT_XRAY] = 1
                    poly.InsertUnder(poly_dict[cryst[i].label])

            # clone to a*b*c grid
            print("clone cell")
            geos = []
            geos_bond = []
            geos_poly = []
            xa, xb, xc = np.meshgrid(np.arange(op[c4d.ID_CIF_A]), np.arange(op[c4d.ID_CIF_B]), np.arange(op[c4d.ID_CIF_C]))
            xa = xa.flatten()
            xb = xb.flatten()
            xc = xc.flatten()

            for i in range(xa.size):
                t = (m[0]*xa[i] + m[1]*xb[i] + m[2]*xc[i])*scale*0.5
                for polynull in rootpoly.GetChildren():
                    # print(polynull.GetName())
                    for p in polynull.GetChildren():
                        newp = p.GetClone(flags=c4d.COPYFLAGS_NONE)
                        newp[c4d.ID_BASEOBJECT_REL_POSITION] = PDTTranslation.nptoc4d(t)
                        geos_poly.append(newp)
                geos.append(PDTFunction.translate(geo,t))
                geos_bond.append(PDTFunction.translate(geo_bond,t))
            
            geo = PDTFunction.merge(geos)
            geo_bond = PDTFunction.merge(geos_bond)
            for poly in geos_poly:
                poly.InsertUnder(poly_dict[poly.GetName()])

            # for label, poly in poly_dict.items():
            #     doc.SetActiveObject(None, c4d.SELECTION_NEW)
            #     print(poly)
            #     children = poly.GetChildren()
            #     for child in children:
            #         print("set active for ", child.GetName())
            #         doc.SetActiveObject(child,c4d.SELECTION_ADD)
            
            # skip replicate site
            print("remove duplicate site")
            arr_atom = np.empty(shape=(0,3))
            for pt in geo.points:
                arr_atom = np.vstack((arr_atom,pt.P))
            arr_atom = arr_atom.astype("int32")
            u, atom_ind = np.unique(arr_atom,axis=0,return_index=True)

            arr_bond = np.empty(shape=(0,3))
            for prim in geo_bond.primitives:
                arr_bond = np.vstack((arr_bond,(geo_bond.points[prim.vertices[0].pointnumber].P + geo_bond.points[prim.vertices[1].pointnumber].P)*0.5))
            arr_bond = arr_bond.astype("int32")
            u, bond_ind = np.unique(arr_bond,axis=0,return_index=True)

            # bond
            sweep = C4DFunction.CreateBond(geo_bond,"label",bond_ind)
            sweep.SetName("bond")
            sweep.InsertUnder(self.rootobj)

            # atom
            rootatom = C4DFunction.CreateAtom(geo,"label",atom_ind)
            rootatom.InsertUnder(self.rootobj)

            # cell
            sweep = C4DFunction.CreateBond(geo_cell,"label")
            sweep.SetName("cell")
            sweep.InsertUnder(self.rootobj)

            self.update = False
            self.hasobj = True
            op.SetDirty(c4d.DIRTYFLAGS_DESCRIPTION)

        if self.hasobj == False:
            return None

        # Change parameters-------------------------------------------------------------
        # read static descriptions
        atomrad = op[c4d.ID_CIF_ATOM_RADIUS_SCALE]
        seg = op[c4d.ID_CIF_ATOM_SEGMENTS]
        bondrad = op[c4d.ID_CIF_BOND_RADIUS]
        showbond = op[c4d.ID_CIF_SHOW_BOND]
        showpoly = op[c4d.ID_CIF_SHOW_POLY]
        showcell = op[c4d.ID_CIF_SHOW_CELL]

        # #change parameters of cloner atom
        clnnull = self.rootobj.GetDown().GetNext().GetDown().GetNext()
        clns = clnnull.GetChildren()
        # clns.reverse()
        for i in range(len(clns)):
            name = self.label_dict[clns[i].GetName()]
            ind = self.symbols.index(name)
            clns[i][c4d.ID_MG_TRANSFORM_COLOR] = self.parameters[2*ind+1]
            clns[i].GetDown()[c4d.PRIM_SPHERE_SUB] = seg
            if op[c4d.ID_CIF_BUILD_MODE] == 1:
                clns[i].GetDown()[c4d.PRIM_SPHERE_RAD] = 10 * bondrad
            else:
                clns[i].GetDown()[c4d.PRIM_SPHERE_RAD] = self.parameters[2*ind]*atomrad
        
        # #change parameters of sweep bond
        sweepnull = self.rootobj.GetDown().GetNext().GetNext()
        sweeps = sweepnull.GetChildren()
        # # print(sweeps)
        # # sweeps.reverse()
        for i in range(len(sweeps)):
            name = self.label_dict[sweeps[i].GetName()]
            ind = self.symbols.index(name)
            sweeps[i].GetDown()[c4d.PRIM_CIRCLE_RADIUS] = 10 * bondrad
            sweeps[i][c4d.ID_BASEOBJECT_COLOR] = self.parameters[2*ind+1]
        # #change parameters of poly
        polynull = self.rootobj.GetDownLast()
        polys = polynull.GetChildren()
        for poly in polys:
            name = self.label_dict[poly.GetName()]
            ind = self.symbols.index(name)
            for p in poly.GetChildren():
                p[c4d.ID_BASEOBJECT_COLOR] = self.parameters[2*ind+1]
        #always build bond and poly, but show it if nessary
        if showcell == 0:
            self.rootobj.GetDown()[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 1   #default is 2, hide is 1, show is 0
            self.rootobj.GetDown()[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 1
        else:
            self.rootobj.GetDown()[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 2   
            self.rootobj.GetDown()[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 2
        
        if showbond == 0:
            sweepnull[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 1   
            sweepnull[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 1
        else:
            sweepnull[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 2
            sweepnull[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 2

        
        if showpoly == 0:
            polynull[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 1
            polynull[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 1
        else:
            polynull[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 2
            polynull[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 2

        return self.rootobj

# main
if __name__ == "__main__":
    # Retrieves the icon path
    directory, _ = os.path.split(__file__)
    fn = os.path.join(directory, "res", "ocifreader.png")

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
                                      g=CHEM4DCrystal,
                                      description="Ocifreader", #file name without .res
                                      info=c4d.OBJECT_GENERATOR|c4d.OBJECT_USECACHECOLOR,
                                      icon=bmp)