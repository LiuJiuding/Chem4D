from typing import Optional
import c4d
from c4d import utils
import scipy.spatial as spt
import pymatgen.core as pmg
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

def FindKValue(center,ligands):
    kt = spt.KDTree(ligands)
    d, ind = kt.query(x=center,k=12)
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
            # print(i, variance)
            if variance>1000.0:
                break
            
    return k-1,kt
        
def ConvexHull(convex_pos:list[c4d.Vector]) -> c4d.PolygonObject:
    hull = spt.ConvexHull(c4dtonp(convex_pos))
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
                        # print("loopstart: ", i)
                        loopstart = i
                if lines[i].find("_atom_site_type_symbol") != -1:
                    type_symbol = i - loopstart -1
                    # print("type_symbol: ", type_symbol)
                if lines[i].find("_atom_site_label") != -1:
                    site_label = i - loopstart - 1
                    # print("site_label: ",site_label)
                if lines[i].find("_atom_site") != -1:
                    if lines[i+1].find("_atom_site") == -1:
                        # print("site_start: ",i+1)
                        site_start = i+1
                        while site_start+site_num < len(lines):
                            if len(lines[site_start+site_num])>5 and lines[site_start+site_num].find("_") == -1:
                                line = lines[site_start+site_num].strip().split()

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
                        # print("site_num: ", site_num)
                        # print("charge_dict: ", charge_dict)
                        # print("label_dict: ", label_dict)
        return cryst,charge_dict,label_dict
 
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
                        # print(str(site.specie))
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

def main() -> None:
    path = "G:\GitClone\Chem4D\sample\LiCoO2.cif"
    a_min = 0
    a_max = 2
    b_min = 0
    b_max = 1
    c_min = 0
    c_max = 1

    cryst, charge_dict, label_dict = ReadCif(path)
    # print(charge_dict)
    # print(cryst.frac_coords)
    # testgeo = PDTGeo()
    # ptnum = PDTFunction.addpoint(testgeo, c4d.Vector(0))
    # PDTFunction.setpointattrib(testgeo,"symbol",0,"Pb")
    # ptnum = PDTFunction.addpoint(testgeo, c4d.Vector(1))
    # # setattr(testgeo.points[0],"new",testgeo.pointattributes["Cd"]())
    # print(testgeo.GetP())
    # print(testgeo.pointattributes["symbol"]())
    anion_label = []
    cation_label = []
    for key,value in charge_dict.items():
        if value == -1:
            anion_label.append(key)
        if value == 1:
            cation_label.append(key)
    
    # set null objects
    rootobj = c4d.BaseObject(5140)
    rootobj[c4d.ID_BASELIST_NAME] = os.path.basename(path)
    rootatom = c4d.BaseObject(5140)
    rootatom[c4d.ID_BASELIST_NAME] = 'atom'
    rootbond = c4d.BaseObject(5140)
    rootbond[c4d.ID_BASELIST_NAME] = 'bond'
    rootpoly = c4d.BaseObject(5140)
    rootpoly[c4d.ID_BASELIST_NAME] = 'poly'
    rootatom.InsertUnder(rootobj)
    rootbond.InsertUnder(rootobj)
    rootpoly.InsertUnder(rootobj)


    # setup pymatgen
    m = cryst.lattice.matrix.T
    print(m)

    cryst.make_supercell([a_max+1,b_max+1,c_max+1],to_unit_cell=True)
    geo = PDTGeo()
    for site in cryst:
        p = site.coords
        f = np.linalg.solve(m,p)
        print(f)
        # a_f = site.frac_coords[0]
        # b_f = site.frac_coords[1]
        # c_f = site.frac_coords[2]
        # pos = site.coords
        # # 剔除目标晶胞范围外的点
        if f[0] <= a_max and f[1] <= b_max and f[2] <= c_max:
            p = nptoc4d(site.coords)*100
            # print(i,j,k,p)
            # print(str(site.specie))
            ptnum = PDTFunction.addpoint(geo, p)
            if type(site.specie)==pmg.periodic_table.Species:
                PDTFunction.setpointattrib(geo, "symbol", ptnum, str(site.specie.element))
            else:
                PDTFunction.setpointattrib(geo, "symbol", ptnum, str(site.specie))
            PDTFunction.setpointattrib(geo, "label", ptnum, site.label)
            PDTFunction.setpointattrib(geo, "charge", ptnum, charge_dict[site.label])
            PDTFunction.setpointgroup(geo, site.label, ptnum, 1)
    labels = geo.pointgroups
    # # 去除所有阴离子

    # geo = BuildPDTGeo(cryst,charge_dict,a_min,a_max,b_min,b_max,c_min,c_max)
    # labels = geo.pointgroups
    # symbols = list(set(geo.GetPointAttrib("symbol")))
    # print(label_dict,anion_label, cation_label)
    # print(labels,symbols)
    # re_anion = geo.GetPointGroupIndex(anion_label)
    # # print(type(re_anion),re_anion)
    # PDTFunction.removepoints(geo, re_anion)
    # # print(labels,symbols,cation_label,anion_label)
    # # 构建一个拓展了阴离子的几何体
    # geo_expand = BuildPDTGeo(cryst,charge_dict,a_min-1,a_max+1,b_min-1,b_max+1,c_min-1,c_max+1)
    # re_cation =  geo_expand.GetPointGroupIndex(cation_label)
    # PDTFunction.removepoints(geo_expand, re_cation)
    # ligands = c4dtonp(geo_expand.GetP())


    # # #收集需要拓展的阴离子位点序号
    # anion_index = []
    # geo_bond = PDT.PDTGeo()
    # for label in cation_label:
    #     centers = c4dtonp(geo.GetP(label))
    #     # 计算cation为中心的凸包多面体和键------------------------------------------------------------------
    #     rootconvex = c4d.BaseObject(5140)
    #     rootconvex[c4d.ID_BASELIST_NAME] = 'poly_'+ label

    #     # calculate k
    #     k,kt = FindKValue(centers[0],ligands)

    #     for center in centers:
    #         # ind是凸包点的序号
    #         d, ind = kt.query(x=center, k=k)
    #         convex_pos = ligands[ind]
    #         anion_index.extend(ind.tolist())
    #         # 生成多面体
    #         if k>3:
    #             poly = ConvexHull(nptoc4d(convex_pos))
    #             poly.SetName(label)
    #             # poly[c4d.ID_BASEOBJECT_USECOLOR] = 2
    #             # poly[c4d.ID_BASEOBJECT_COLOR] = read_dprop[self.label_dict[label]][1]
    #             # poly[c4d.ID_BASEOBJECT_XRAY] = 1
    #             poly.InsertUnder(rootconvex)

    #         # 构筑化学键
    #         for i in ind:
    #             beginpos = nptoc4d(center)
    #             endpos = nptoc4d(ligands[i])
    #             midpos = (beginpos + endpos) * 0.5

    #             beginnum = PDTFunction.addpoint(geo_bond,beginpos)
    #             midnum = PDTFunction.addpoint(geo_bond,midpos)
    #             endnum = PDTFunction.addpoint(geo_bond,endpos)

    #             bondbegin = PDTFunction.addprim(geo_bond,'polyline',[beginnum,midnum])
    #             bondend = PDTFunction.addprim(geo_bond,'polyline',[midnum,endnum])

    #             beginsymbol = label_dict[label]
    #             endsymbol = getattr(geo_expand.points[i],"symbol")
    #             PDTFunction.setprimattrib(geo_bond,'symbol',bondbegin,beginsymbol)
    #             PDTFunction.setprimattrib(geo_bond,'symbol',bondend,endsymbol)
    #     rootconvex.InsertUnder(rootpoly)
    # # print(len(geo_bond.primitives[0].vertices))
    # # print(labels,symbols)
    # for symbol in symbols:
    #     sweep = C4DFunction.CreateSweep(geo_bond,symbol)
    #     # sweep[c4d.ID_BASEOBJECT_USECOLOR] = 2
    #     # sweep[c4d.ID_BASEOBJECT_COLOR] = self.read_dprop[symbol][1]
    #     sweep.InsertUnder(rootbond)

    # # Generate atoms-------------------------------------------------------------------------------
    # # Merge anion to geo
    # for i in anion_index:
    #     ptnum = PDTFunction.addpoint(geo, nptoc4d(ligands[i]))
    #     PDTFunction.setpointattrib(geo, "symbol", ptnum, getattr(geo_expand.points[i],"symbol"))
    #     PDTFunction.setpointattrib(geo, "label", ptnum, getattr(geo_expand.points[i],"label"))
    #     PDTFunction.setpointattrib(geo, "charge", ptnum, getattr(geo_expand.points[i],"charge"))
    #     PDTFunction.setpointgroup(geo, getattr(geo_expand.points[i],"label"), ptnum, 1)
    # Create mesh
    mesh = C4DFunction.CreateMesh(geo)
    mesh.InsertUnder(rootobj)
    # Create sph and cloner
    for label in labels:
        sph = C4DFunction.AddSphere(label,50)
        # print(nptoc4d(np.random.randn(1,3)))
        cln = C4DFunction.AddCloner(mesh, label, label, nptoc4d(np.random.randn(1,3))[0])
        sph.InsertUnder(cln)
        cln.InsertUnder(rootatom)

    doc.InsertObject(rootobj)
    c4d.EventAdd()
if __name__ == '__main__':
    main()
