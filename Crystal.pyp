################# Chem4D Crystal V1.0.0 ##################
from typing import Optional
import c4d
from c4d import utils
import os
import copy

import gemmi
import pymatgen.core as pmg

import numpy as np
import scipy.spatial as spt

import PDT
from PDT import PDTGeo, PDTFunction, C4DFunction

PLUGINID = 	1061833
PLUGINNAME = "Chem4D Crystal"
VERSION = "1.0.0"
TITLE = f"{PLUGINNAME} v{VERSION}"
HELP = "crystal building tool"
# Dynamic group and parameters IDs
CHEM4D_DYNAMICGROUP = 1100
CHEM4D_DYNAMICGROUP_FIRSTPARAMETER = CHEM4D_DYNAMICGROUP + 1


def nptoc4d(nparray:np.array):
    if nparray.ndim == 1:
        return c4d.Vector(nparray[0],nparray[1],nparray[2])
    else:
        c4dlist = []
        for i in nparray:
            c4dlist.append(c4d.Vector(i[0],i[1],i[2]))
        return c4dlist


def c4dtonp(c4darray):
    if type(c4darray)==c4d.Vector:
        return np.array([c4darray[0],c4darray[1],c4darray[2]])
    else:
        newlist = []
        for i in c4darray:
            newlist.append([i[0],i[1],i[2]])
        return np.array(newlist)
        
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
        # print("read property ok")
        return dprop
    
    @staticmethod
    def FindKValue(center,ligands):
        kt = spt.KDTree(ligands)
        d, ind = kt.query(x=center,k=9)
        distance = []
        for i in ind:
            distance.append(np.linalg.norm(ligands[i]-center))
        distance.sort()
        # print(distance)
        k:int
        for i in range(len(distance)):
            if i>3:
                sample = distance[:i]
                variance = np.var(sample)
                k = i
                print(i, variance)
                if variance>1000.0:
                    break
                
        return k-1,kt
    @staticmethod
    def ConvexHull(convex_pos:list[c4d.Vector]) -> c4d.PolygonObject:
        hull = spt.ConvexHull(c4dtonp(convex_pos))
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
    
    @staticmethod
    def ReadCif(path):
        if path == '':
            return False
        if os.path.isfile(path) == False:
            print('path is not valid')
            return False
        #get crystal
        cryst = pmg.Structure.from_file(path)

        if cryst is None or False:
            print('failed reading .cif file')
            return False
        else:
            print('successfully reading cif file\n',path)

        with open(path,'r') as f:
            lines = f.readlines()
            loopstart:int = 0
            type_symbol:int = 0
            site_label:int = 0
            site_start:int = 0
            site_num:int = 0
            charge_dict = {}
            label_dict = {}
            for i in range(len(lines)):
                if lines[i].find("loop_") != -1:
                    if lines[i+1].find("_atom_site") != -1:
                        print("loopstart: ", i)
                        loopstart = i
                if lines[i].find("_atom_site_type_symbol") != -1:
                    type_symbol = i - loopstart -1
                    print("type_symbol: ", type_symbol)
                if lines[i].find("_atom_site_label") != -1:
                    site_label = i - loopstart - 1
                    print("site_label: ",site_label)
                if lines[i].find("_atom_site") != -1:
                    if lines[i+1].find("_atom_site") == -1:
                        print("site_start: ",i+1)
                        site_start = i+1
                        while site_start+site_num < len(lines):
                            if len(lines[site_start+site_num])>5 and lines[site_start+site_num].find("_") == -1:
                                line = lines[site_start+site_num].strip().split()
                                print(line)
                                if len(line[type_symbol]) in {1,2}:
                                    charge_dict[line[site_label]] = 0
                                    label_dict[line[site_label]] = line[site_label]
                                elif line[type_symbol][-2] == "+" or line[type_symbol][-1] == "+":
                                    charge_dict[line[site_label]] = 1
                                    label_dict[line[site_label]] = line[type_symbol][:-2]
                                elif line[type_symbol][-2] == "-" or line[type_symbol][-1] == "-":
                                    charge_dict[line[site_label]] = -1
                                    label_dict[line[site_label]] = line[type_symbol][:-2]
                            site_num = site_num + 1
                        print("site_num: ", site_num)
                        print("charge_dict: ", charge_dict)
                        print("label_dict: ", label_dict)
        return cryst,charge_dict,label_dict
        
    @staticmethod
    def BuildPDTGeo(cryst, charge_dict, a_min=0, a_max=1, b_min=0, b_max=1, c_min=0, c_max=1):
        a = cryst.lattice.matrix[0]
        b = cryst.lattice.matrix[1]
        c = cryst.lattice.matrix[2]
        geo = PDTGeo()
        for site in cryst:
            a_f = site.frac_coords[0]
            b_f = site.frac_coords[1]
            c_f = site.frac_coords[2]
            pos = site.coords
            # cryst中的原子是经过对称变换的，平移后不会有重复点，只需扩展后剔除目标晶胞范围外的点
            for i in range(a_min,a_max+1):
                for j in range(b_min,b_max+1):
                    for k in range(c_min,c_max+1):
                        if a_f+i<=a_max and b_f+j<=b_max and c_f+k<=c_max:
                            p = nptoc4d(pos+a*i + b*j + c*k)*100
                            # print(i,j,k,p)
                            ptnum = PDTFunction.addpoint(geo, p)
                            if type(site.specie)==pmg.periodic_table.Species:
                                PDTFunction.setpointattrib(geo, "symbol", ptnum, str(site.specie.element))
                            else:
                                PDTFunction.setpointattrib(geo, "symbol", ptnum, str(site.specie))
                            PDTFunction.setpointattrib(geo, "label", ptnum, site.label)
                            PDTFunction.setpointattrib(geo, "charge", ptnum, charge_dict[site.label])
                            PDTFunction.setpointgroup(geo, site.label, ptnum, 1)
                            # PDTFunction.setpointgroup(geo, atm.element, ptnum, 1)
        return geo

class CHEM4DCrystal(c4d.plugins.ObjectData, CHEM4DHelper):

    def __init__(self):
        self.hasobj = False
        self.update = False
        self.read = False
        self.read_dprop = {}
        self.parameters = []
        self.geo = PDTGeo()

        self.symbols = []
        self.labels = []
        self.label_dict = {}
        # self.geo_bond = PDTGeo()
        self.rootobj = c4d.BaseObject(5140)
        self.SetOptimizeCache(True)

    def Init(self, node):
        node[c4d.ID_CIF_ATOM_RADIUS_SCALE] = 1.0
        node[c4d.ID_CIF_ATOM_SEGMENTS] = 20
        node[c4d.ID_CIF_BOND_RADIUS] = 1.0
        node[c4d.ID_CIF_BUILD_MODE] = 0
        node[c4d.ID_CIF_A_MIN] = 0
        node[c4d.ID_CIF_A_MAX] = 1
        node[c4d.ID_CIF_B_MIN] = 0
        node[c4d.ID_CIF_B_MAX] = 1
        node[c4d.ID_CIF_C_MIN] = 0
        node[c4d.ID_CIF_C_MAX] = 1

        #read property file
        directory, _ = os.path.split(__file__)
        fn = os.path.join(directory, "elements.txt")
        self.read_dprop = self.ReadProperty(fn)
        return True

    def Read(self, node, hf, level):
        # Reads the number of dynamic parameters
        node[c4d.ID_CIF_ATOM_RADIUS_SCALE] = hf.ReadFloat32()
        node[c4d.ID_CIF_ATOM_SEGMENTS] = hf.ReadInt32()
        node[c4d.ID_CIF_BOND_RADIUS] = hf.ReadFloat32()
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
        # dest.geo = copy.copy(self.geo)
        # dest.geo_bond = copy.copy(self.geo_bond)
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
                self.parameters.append(self.read_dprop[i][0])
                self.parameters.append(self.read_dprop[i][1])

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
        if self.update == True:
            self.rootobj = c4d.BaseObject(5140)
            a_min = op[c4d.ID_CIF_A_MIN]
            a_max = op[c4d.ID_CIF_A_MAX]
            b_min = op[c4d.ID_CIF_B_MIN]
            b_max = op[c4d.ID_CIF_B_MAX]
            c_min = op[c4d.ID_CIF_C_MIN]
            c_max = op[c4d.ID_CIF_C_MAX]
            # check path, cif file and creat crystals cryst-------------------------------------------------------
            cryst, charge_dict, label_dict = self.ReadCif(op[c4d.ID_CIF_PATH])
            anion_label = []
            cation_label = []
            for key,value in charge_dict.items():
                if value == -1:
                    anion_label.append(key)
                if value == 1:
                    cation_label.append(key)
            if cryst == False:
                self.hasobj = False
                self.update = False
                return self.rootobj
            if self.read == False:
                self.parameters.clear()
                for i in self.symbols:
                    self.parameters.append(self.read_dprop[i][0])
                    self.parameters.append(self.read_dprop[i][1])
                self.update = False
                self.hasobj = True
            if a_min<a_max and b_min<b_max and c_min<c_max:
                # set null objects------------------------------------------------------------------------------
                self.rootobj[c4d.ID_BASELIST_NAME] = os.path.basename(op[c4d.ID_CIF_PATH])
                rootatom = c4d.BaseObject(5140)
                rootatom[c4d.ID_BASELIST_NAME] = 'atom'
                rootbond = c4d.BaseObject(5140)
                rootbond[c4d.ID_BASELIST_NAME] = 'bond'
                rootpoly = c4d.BaseObject(5140)
                rootpoly[c4d.ID_BASELIST_NAME] = 'poly'
                rootatom.InsertUnder(self.rootobj)
                rootbond.InsertUnder(self.rootobj)
                rootpoly.InsertUnder(self.rootobj)
                # 构建目标大小的几何体------------------------------------------------------------
                self.geo = self.BuildPDTGeo(cryst,charge_dict,a_min,a_max,b_min,b_max,c_min,c_max)
                self.labels = self.geo.pointgroups
                self.symbols = list(set(self.geo.GetPointAttrib("symbol")))
                self.label_dict = label_dict
                PDTFunction.removepoints(self.geo, self.geo.GetPointGroupIndex(anion_label))
                
                # 构建一个拓展了阴离子的几何体
                geo_expand = self.BuildPDTGeo(cryst,charge_dict,a_min-1,a_max+1,b_min-1,b_max+1,c_min-1,c_max+1)
                PDTFunction.removepoints(geo_expand, geo_expand.GetPointGroupIndex(cation_label))
                ligands = c4dtonp(geo_expand.GetP())

                # 收集需要拓展的阴离子位点序号
                anion_index = []
                geo_bond = PDT.PDTGeo()
                for label in cation_label:
                    centers = c4dtonp(self.geo.GetP(label))

                    # 计算cation为中心的凸包多面体和键------------------------------------------------------------------
                    rootconvex = c4d.BaseObject(5140)
                    rootconvex[c4d.ID_BASELIST_NAME] = 'poly_'+ label

                    # calculate k
                    k,kt = self.FindKValue(centers[0],ligands)
                    
                    for center in centers:
                        # ind是凸包点的序号
                        d, ind = kt.query(x=center, k=k)
                        convex_pos = ligands[ind]
                        anion_index.extend(ind.tolist())
                        # 生成多面体
                        if k>3:
                            poly = self.ConvexHull(nptoc4d(convex_pos))
                            poly.SetName(label)
                            poly[c4d.ID_BASEOBJECT_USECOLOR] = 2
                            poly[c4d.ID_BASEOBJECT_COLOR] = self.read_dprop[self.label_dict[label]][1]
                            poly[c4d.ID_BASEOBJECT_XRAY] = 1
                            poly.InsertUnder(rootconvex)

                        # 构筑化学键
                        for i in ind:
                            beginpos = nptoc4d(center)
                            endpos = nptoc4d(ligands[i])
                            midpos = (beginpos + endpos) * 0.5

                            beginnum = PDTFunction.addpoint(geo_bond,beginpos)
                            midnum = PDTFunction.addpoint(geo_bond,midpos)
                            endnum = PDTFunction.addpoint(geo_bond,endpos)

                            bondbegin = PDTFunction.addprim(geo_bond,'polyline',[beginnum,midnum])
                            bondend = PDTFunction.addprim(geo_bond,'polyline',[midnum,endnum])

                            beginsymbol = self.label_dict[label]
                            endsymbol = getattr(geo_expand.points[i],"symbol")
                            PDTFunction.setprimattrib(geo_bond,'symbol',bondbegin,beginsymbol)
                            PDTFunction.setprimattrib(geo_bond,'symbol',bondend,endsymbol)


                    rootconvex.InsertUnder(rootpoly)
                for symbol in self.symbols:
                    sweep = C4DFunction.CreateSweep(geo_bond,symbol)
                    sweep[c4d.ID_BASEOBJECT_USECOLOR] = 2
                    sweep[c4d.ID_BASEOBJECT_COLOR] = self.read_dprop[symbol][1]
                    sweep.InsertUnder(rootbond)

                # Generate atoms-------------------------------------------------------------------------------
                # Create mesh
                for i in anion_index:
                    ptnum = PDTFunction.addpoint(self.geo, nptoc4d(ligands[i]))
                    PDTFunction.setpointattrib(self.geo, "symbol", ptnum, getattr(geo_expand.points[i],"symbol"))
                    PDTFunction.setpointattrib(self.geo, "label", ptnum, getattr(geo_expand.points[i],"label"))
                    PDTFunction.setpointattrib(self.geo, "charge", ptnum, getattr(geo_expand.points[i],"charge"))
                    PDTFunction.setpointgroup(self.geo, getattr(geo_expand.points[i],"label"), ptnum, 1)
                mesh = C4DFunction.CreateMesh(self.geo)
                mesh.InsertUnder(self.rootobj)
                # Create sph and cloner
                for label in self.labels:
                    sph = C4DFunction.AddSphere(label, self.read_dprop[self.label_dict[label]][0])
                    cln = C4DFunction.AddCloner(mesh, label, label, self.read_dprop[self.label_dict[label]][1])
                    sph.InsertUnder(cln)
                    cln.InsertUnder(rootatom)

                self.update = False
                self.hasobj = True

        if self.hasobj == False:
                return None
        if self.hasobj == True:   
            #Change parameters-------------------------------------------------------------
            #read static descriptions
            # atomrad = op[c4d.ID_CIF_ATOM_RADIUS_SCALE]
            # seg = op[c4d.ID_CIF_ATOM_SEGMENTS]
            # bondrad = op[c4d.ID_CIF_BOND_RADIUS]


            # #change parameters of cloner atom
            # clns = self.rootobj.GetDownLast().GetChildren()
            # # clns.reverse()
            # for i in range(len(clns)):
            #     name = self.realsymbol[clns[i].GetName()]
            #     ind = self.symbols.index(name)
            #     clns[i][c4d.ID_MG_TRANSFORM_COLOR] = self.parameters[2*ind+1]
            #     clns[i].GetDown()[c4d.PRIM_SPHERE_SUB] = seg
            #     if op[c4d.ID_CIF_BUILD_MODE] == 3:
            #         clns[i].GetDown()[c4d.PRIM_SPHERE_RAD] = 10 * bondrad
            #     else:
            #         clns[i].GetDown()[c4d.PRIM_SPHERE_RAD] = self.parameters[2*ind]*atomrad
            
            # #change parameters of sweep bond
            # sweeps = self.rootobj.GetDown().GetNext().GetNext().GetChildren()
            # # print(sweeps)
            # # sweeps.reverse()
            # for i in range(len(sweeps)):
            #     name = sweeps[i].GetName()
            #     ind = self.symbols.index(name)
            #     sweeps[i].GetDown()[c4d.PRIM_CIRCLE_RADIUS] = 10 * bondrad
            #     sweeps[i][c4d.ID_BASEOBJECT_COLOR] = self.parameters[2*ind+1]
            
            # #always build bond, but show it if nessary
            # if op[c4d.ID_MOL_BUILD_MODE] == 1:
            #     self.rootobj.GetDown().GetNext()[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 1   #default is 2, hide is 1, show is 0
            #     self.rootobj.GetDown().GetNext()[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 1
            # else:
            #     self.rootobj.GetDown().GetNext()[c4d.ID_BASEOBJECT_VISIBILITY_EDITOR] = 2
            #     self.rootobj.GetDown().GetNext()[c4d.ID_BASEOBJECT_VISIBILITY_RENDER] = 2
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