from typing import Optional
import c4d
from c4d import utils
import scipy.spatial as spt
import pymatgen.core as pmg
import pymatgen.analysis.local_env as pmg_env
import numpy as np

import PDT
import os
from PDT import PDTGeo, PDTFunction, C4DFunction

doc: c4d.documents.BaseDocument  # The active document
op: Optional[c4d.BaseObject]  # The active object, None if unselected


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

def getsymbol(site):
    symbol = ""
    if type(site.specie)==pmg.periodic_table.Species:
        symbol = str(site.specie.element)
    else:
        symbol = str(site.specie)
    return symbol

def ReadProperty(file):
    fprop = open(file,'r')
    dprop = {}
    for line in fprop:
        f2 = line.strip().split()
        dprop[f2[0]] = [int(f2[1]),c4d.Vector(int(f2[2]),int(f2[3]),int(f2[4]))/255]   #如果1个key有多个value，可以写成new_dict[f2[0]]  =  f2[1:]

    fprop.close()
    # print("read property ok")
    return dprop   

def ConvexHull(convex_pos:list[c4d.Vector]) -> c4d.PolygonObject:
    makeconvex = True
    try:
        hull = spt.ConvexHull(c4dtonp(convex_pos))
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

def BuildCrystal(path,dprop):

    scale = 100
    delta = 0.0001
    a = 1
    b = 1
    c = 1
    # read cif, get dictionary and label list----------------------------------------
    cryst = ReadCif(path)
    cryst.add_oxidation_state_by_guess()
    charge_dict = {}
    label_dict = {}
    for site in cryst:
        charge_dict[site.label] = str(site.specie)[-1]
        label_dict[site.label] = str(site.specie.element)

    anion_label = []
    cation_label = []
    for key,value in charge_dict.items():
        if value == "-":
            anion_label.append(key)
        if value == "+":
            cation_label.append(key)
    # print(anion_label,cation_label)
    # set null objects----------------------------------------------------------------
    rootobj = c4d.BaseObject(5140)
    rootobj[c4d.ID_BASELIST_NAME] = os.path.basename(path)
    # rootatom = c4d.BaseObject(5140)
    # rootatom[c4d.ID_BASELIST_NAME] = 'atom'
    # rootbond = c4d.BaseObject(5140)
    # rootbond[c4d.ID_BASELIST_NAME] = 'bond'
    rootpoly = c4d.BaseObject(5140)
    rootpoly[c4d.ID_BASELIST_NAME] = 'poly'
    
    # rootbond.InsertUnder(rootobj)
    rootpoly.InsertUnder(rootobj)

    # setup unit cell line-------------------------------------------------------------------
    geo_cell = PDTGeo()
    axises = ["a","b","c"]
    colors = [c4d.Vector(1,0,0),c4d.Vector(0,1,0),c4d.Vector(0,0,1)]
    for i in range(3):
        beginpos = c4d.Vector(0)
        endpos = nptoc4d(cryst.lattice.matrix[i])*scale
        beginnum = PDTFunction.addpoint(geo_cell,beginpos)
        endnum = PDTFunction.addpoint(geo_cell,endpos)
        celledge = PDTFunction.addprim(geo_cell,'polyline',[beginnum,endnum])
        PDTFunction.setprimattrib(geo_cell,'label',celledge,axises[i])
        PDTFunction.setprimattrib(geo_cell,'Cd',celledge,colors[i])

    # make super cell to fill unit cell----------------------------------------------------
    cryst.make_supercell([a+1,b+1,c+1])
    m = cryst.lattice.matrix

    geo = PDTGeo()
    NN = pmg_env.CrystalNN(distance_cutoffs=None,x_diff_weight=0,porous_adjustment=False)
    # build atom inside boundary-----------------------------------------------------------
    inside = []
    for i in range(len(cryst)):

        minus_a = 0
        minus_b = 0
        minus_c = 0
        f = cryst[i].frac_coords

        if f[0] <= 1-1/(a+1)+delta and f[1] <= 1-1/(b+1)+delta and f[2] <= 1-1/(c+1)+delta:

            inside.append(i)
            p = nptoc4d(cryst[i].coords)*scale
            ptnum = PDTFunction.addpoint(geo, p)
            
            symbol = getsymbol(cryst[i])

            PDTFunction.setpointattrib(geo, "symbol", ptnum, symbol)
            PDTFunction.setpointattrib(geo, "label", ptnum, cryst[i].label)
            PDTFunction.setpointattrib(geo, "Cd", ptnum, dprop[symbol][1])
            PDTFunction.setpointattrib(geo, "charge", ptnum, charge_dict[cryst[i].label])
            PDTFunction.setpointgroup(geo, cryst[i].label, ptnum, 1)
    labels = geo.pointgroups

    # build poly and bond-------------------------------------------------------------------
    geo_bond = PDTGeo()
    poly_dict = {}
    for label in cation_label:
        rootconvex = c4d.BaseObject(5140)
        rootconvex[c4d.ID_BASELIST_NAME] = 'poly_'+ label
        rootconvex.InsertUnder(rootpoly)
        poly_dict[label] = rootconvex
    for i in inside:
        dicts = NN.get_nn_info(cryst,i)
        # begin property
        beginpos = nptoc4d(cryst[i].coords)*scale
        beginlabel = cryst[i].label
        beginsymbol = getsymbol(cryst[i])
        pos = []
        for dic in dicts:
            # end property
            endpos = nptoc4d(dic["site"].coords)*scale
            pos.append(endpos)
            endlabel = dic["site"].label
            endsymbol = getsymbol(dic["site"])

            midpos = (beginpos + endpos) * 0.5

            p = dic["site"].coords
            f = np.linalg.solve(m.T,p)
            
            outside = f[0] > float(1-1/(a+1)+delta) or f[1] > float(1-1/(b+1)+delta) or f[2] > float(1-1/(c+1)+delta) or f[0] < 0.0 or f[1] < 0.0 or f[2] < 0.0
            # print(f,outside)
            # 创建中心bond：中心是阴离子 且 配体在边界之外，不创建中心bond
            if not(charge_dict[cryst[i].label] == "-" and outside):
                # print(cryst[i].label,f)
                beginnum = PDTFunction.addpoint(geo_bond,beginpos)
                midnum = PDTFunction.addpoint(geo_bond,midpos)
                bondbegin = PDTFunction.addprim(geo_bond,'polyline',[beginnum,midnum])
                PDTFunction.setprimattrib(geo_bond,'label',bondbegin,beginlabel)
                PDTFunction.setprimattrib(geo_bond,'Cd',bondbegin,dprop[beginsymbol][1])

            # 创建配体atom和bond：中心是阳离子 且 配体在边界之外，创建配体atom和配体bond
            if charge_dict[cryst[i].label] == "+" and outside:
                # atom
                # print(f,endpos)
                ptnum = PDTFunction.addpoint(geo, endpos)
                PDTFunction.setpointattrib(geo, "symbol", ptnum, endsymbol)
                PDTFunction.setpointattrib(geo, "label", ptnum, endlabel)
                PDTFunction.setpointattrib(geo, "Cd", ptnum, dprop[endsymbol][1])
                PDTFunction.setpointattrib(geo, "charge", ptnum, charge_dict[endlabel])
                PDTFunction.setpointgroup(geo, endlabel, ptnum, 1)
                # bond
                endnum = PDTFunction.addpoint(geo_bond,endpos)
                midnum = PDTFunction.addpoint(geo_bond,midpos)
                bondend = PDTFunction.addprim(geo_bond,'polyline',[endnum,midnum])
                PDTFunction.setprimattrib(geo_bond,'label',bondend,endlabel)
                PDTFunction.setprimattrib(geo_bond,'Cd',bondend,dprop[endsymbol][1])

        # print(len(pos))
        if charge_dict[cryst[i].label] == "+" and len(pos)>3:
            poly = ConvexHull(pos)
            poly.SetName(cryst[i].label)
            poly[c4d.ID_BASEOBJECT_USECOLOR] = 2
            poly[c4d.ID_BASEOBJECT_COLOR] = dprop[getsymbol(cryst[i])][1]
            poly[c4d.ID_BASEOBJECT_XRAY] = 1
            poly.InsertUnder(poly_dict[cryst[i].label])
    # print(geo)
    # print(geo_bond)
    sweep = C4DFunction.CreateBond(geo_bond,"label")
    sweep.SetName("bond")
    # sweep[c4d.ID_BASEOBJECT_USECOLOR] = 2
    # sweep[c4d.ID_BASEOBJECT_COLOR] = self.read_dprop[symbol][1]
    sweep.InsertUnder(rootobj)
    # cell 
    sweep = C4DFunction.CreateBond(geo_cell,"label")
    sweep.SetName("cell")
    # sweep[c4d.ID_BASEOBJECT_USECOLOR] = 2
    # sweep[c4d.ID_BASEOBJECT_COLOR] = self.read_dprop[symbol][1]
    sweep.InsertUnder(rootobj)

    # Create mesh
    rootatom = C4DFunction.CreateAtom(geo)
    rootatom.InsertUnder(rootobj)


    doc.InsertObject(rootobj)
    c4d.EventAdd()

def main() -> None:
    # path = "G:\GitClone\Chem4D\sample\Fe(CO)5.cif"
    path = "G:\GitClone\Chem4D\sample\LiCoO2.cif"
    
    #read property file
    directory, _ = os.path.split(__file__)
    fn = os.path.join(r"G:\GitClone\Chem4D", "elements.txt")
    dprop = ReadProperty(fn)
    BuildCrystal(path,dprop)
    # targetdir = r"G:\GitClone\Chem4D\sample\test_structures\common_binaries"
    # files = os.listdir(targetdir)
    # for file in files:
    #     BuildCrystal(os.path.join(targetdir,file),dprop)
    
if __name__ == '__main__':
    main()
