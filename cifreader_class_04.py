from typing import Optional
import c4d
from c4d import utils
import scipy.spatial as spt
from crystals import Crystal
import numpy as np
import PDT
import gemmi
from PDT import PDTFunction, C4DFunction
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

def ReadCif(path, a_min=0, a_max=1, b_min=0, b_max=1, c_min=0, c_max=1):
    cryst = Crystal.from_cif(path)
    groups = cryst.groupby(by="equivalent_atoms")
    print(type(groups))

    print(groups)
    a,b,c = cryst.lattice_vectors

    geo = PDT.PDTGeo()
    for atm in cryst:
        a_f = atm.coords_fractional[0]
        b_f = atm.coords_fractional[1]
        c_f = atm.coords_fractional[2]

        symbol = atm.element
        pos = nptoc4d(atm.coords_cartesian)
        # cryst中的原子是经过对称变换的，平移后不会有重复点，只需扩展后剔除目标晶胞范围外的点
        for i in range(a_min,a_max+1):
            for j in range(b_min,b_max+1):
                for k in range(c_min,c_max+1):
                    if a_f+i<=a_max and b_f+j<=b_max and c_f+k<=c_max:
                        ptnum = PDTFunction.addpoint(geo, (pos + nptoc4d(a*i + b*j + c*k)) *100)
                        PDTFunction.setpointattrib(geo, "symbol", ptnum, symbol)
                        PDTFunction.setpointgroup(geo, symbol, ptnum, 1)
    return geo

def main() -> None:
    path = "G:\GitClone\Chem4D\sample\Li2MnO3.cif"
    a_min = 0
    a_max = 1
    b_min = 0
    b_max = 1
    c_min = 0
    c_max = 1

    cif_doc = gemmi.cif.read(path)
    LCO = gemmi.make_small_structure_from_block(cif_doc.sole_block())

    geo = ReadCif(path,a_min,a_max,b_min,b_max,c_min,c_max)

    cations = []
    anions = []
    isatom = False
    for site in LCO.sites:
        # print(site.type_symbol)
        if site.type_symbol[-1] == "+":
            cations.append(site.type_symbol[:-2])
        elif site.type_symbol[-1] == "-": 
            anions.append(site.type_symbol[:-2])
        elif site.type_symbol[-1] == "0": 
            cations.append(site.type_symbol[:-1])
            anions.append(site.type_symbol[:-1])
            isatom = True

    # print(cations,anions)
    
    # 去除所有阴离子
    # if not isatom:
    for anion in anions:
        idx = geo.GetPointGroupIndex(anion)
        PDTFunction.removepoints(geo, idx)
    # 重新生成拓展的阴离子
    geo_expand = ReadCif(path,a_min-1,a_max+1,b_min-1,b_max+1,c_min-1,c_max+1)

    rootobj = c4d.BaseObject(5140)
    rootobj[c4d.ID_BASELIST_NAME] = 'rootobj'
    #收集需要拓展的阴离子位点序号
    anion_index = []
    for cation in cations:
        centers = c4dtonp(geo.GetP(cation))
        ligands = c4dtonp(geo_expand.GetP(anions[0]))



        # 计算cation为中心的凸包多面体和键------------------------------------------------------------------
        rootpoly = c4d.BaseObject(5140)
        rootpoly[c4d.ID_BASELIST_NAME] = 'poly'+cation
        kt = spt.KDTree(ligands)

        # 凸包的k取值算法
        # 1.取1*1和3*3晶胞里的阳离子和阴离子分别做中心和配体（减小计算量）or 直接计算k=12的凸包点
        # 2.计算12个距离并从小到大排序
        # 3.用i+1的减i距离，若距离差大于threshold值，则k=i+1
        d, ind = kt.query(x=centers[1],k=8)
        distance = []
        for i in ind:
            distance.append(np.linalg.norm(ligands[i]-centers[1]))
        distance.sort()
        # print(distance)
        # if isatom:
        #     distance.pop[0]
        k: int
        for i in range(len(distance)-1):
            # print(i)
            if distance[i+1]-distance[i]>100:
                k = i+1
                break
        # print(k)
        geo_bond = PDT.PDTGeo()
        for center in centers:
            # ind是凸包点的序号
            d, ind = kt.query(x=center, k=k)
            convex_pos = ligands[ind]
            # print(ind)
            anion_index.extend(ind.tolist())
            # 生成多面体
            if k>3:
                poly = ConvexHull(nptoc4d(convex_pos))
                poly.InsertUnder(rootpoly)

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

                beginsymbol = cation
                endsymbol = anion[0]
                PDTFunction.setprimattrib(geo_bond,'symbol',bondbegin,beginsymbol)
                PDTFunction.setprimattrib(geo_bond,'symbol',bondend,endsymbol)

        rootpoly.InsertUnder(rootobj)
        sweepO = C4DFunction.CreateSweep(geo_bond,anion[0])
        sweepO.InsertUnder(rootobj)
        sweepP = C4DFunction.CreateSweep(geo_bond,cation)
        sweepP.InsertUnder(rootobj)


    # Generate atoms
    # print(anion_index)
    anion_index = np.unique(anion_index)
    for i in anion_index:
        ptnum = PDTFunction.addpoint(geo, nptoc4d(ligands[i]))
        PDTFunction.setpointattrib(geo, "symbol", ptnum, anions[0])
        PDTFunction.setpointgroup(geo, anions[0], ptnum, 1)

    mesh = C4DFunction.CreateMesh(geo)
    mesh.InsertUnder(rootobj)
    for element in geo.pointgroups:
        sph = C4DFunction.AddSphere(element, 50)
        cln = C4DFunction.AddCloner(mesh, element, element)
        sph.InsertUnder(cln)
        cln.InsertUnder(rootobj)
    doc.InsertObject(rootobj)
    c4d.EventAdd()

if __name__ == '__main__':
    main()
