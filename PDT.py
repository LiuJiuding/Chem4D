import c4d
import copy
import numpy as np


class PDTTranslation():
    @staticmethod
    def listtoc4d(mylist:list):
        c4dlist = []
        for i in mylist:
            c4dlist.append(c4d.Vector(i[0],i[1],i[2]))
        return c4dlist
        
    @staticmethod
    def nptoc4d(nparray:np.ndarray):
        if nparray.ndim == 1:
            return c4d.Vector(nparray[0],nparray[1],nparray[2])
        else:
            c4dlist = []
            for i in nparray:
                c4dlist.append(c4d.Vector(i[0],i[1],i[2]))
            return c4dlist
    @staticmethod
    def c4dtonp(c4darray):
        if type(c4darray)==c4d.Vector:
            return np.array([c4darray[0],c4darray[1],c4darray[2]])
        else:
            newlist = []
            for i in c4darray:
                newlist.append([i[0],i[1],i[2]])
            return np.array(newlist)

# attribute type: int, str, float, np.ndarray(3,)
# point and primitive attribute dictionary stores type for int, str and float, stores string for ndarray
class PDTPoint(object):
    def __init__(self, pos:np.ndarray = np.array([0.0,0.0,0.0])) -> None:
        self.P = pos
class PDTVertice(object):
    def __init__(self, ptnum:int) -> None:
        self.pointnumber = ptnum
class PDTPrimitive(object):
    def __init__(self,ptnums:list, primtype:str) -> None:
        self.primtype = primtype
        self.pointcount = len(ptnums)
        self.vertices = []
        for i in ptnums:
            v = PDTVertice(i)
            self.vertices.append(v)
class PDTGeo(object):
    def __init__(self, pts:np.ndarray = np.empty(shape=(0,3))) -> None:
        self.points = []
        for i in pts:
            self.points.append(PDTPoint(i))   
        self.primitives = []
        self.vertices = []

        self.pointcount = len(pts)
        self.primitivecount = len(self.primitives)

        self.pointattributes = {"P":"ndarray3"}
        self.primitiveattributes = {}

        self.pointgroups = []
    # 提取点位置列表，可以按组名提取
    def GetP(self,groupname='') -> list:
        pos_list = []
        if len(groupname) == 0:
            for pt in self.points:
                pos_list.append(pt.P)
            return PDTTranslation.listtoc4d(pos_list)
        
        elif type(groupname) == str:
            if groupname in self.pointgroups:
                for pt in self.points:
                    if getattr(pt,'group:'+groupname) == 1:
                        pos_list.append(pt.P)
                return PDTTranslation.listtoc4d(pos_list)
            else:
                print('group '+ groupname + ' does not exist')
                return False
        elif type(groupname) == list:
            for name in groupname:
                if name in self.pointgroups:
                    for pt in self.points:
                        if getattr(pt,'group:'+name) == 1:
                            pos_list.append(c4d.Vector(pt.P))
                else:
                    print('group '+ name + ' does not exist')
                    return False
            return PDTTranslation.listtoc4d(pos_list)
        else:
            print("invalid input name")
            return False
    # 提取点属性列表
    def GetPointAttrib(self, name:str, index:int = -1):
        if name in self.pointattributes:
            if index in range(self.pointcount):
                # print("single index")
                attrib = getattr(self.points[index],name)
                return attrib
            else:
                # print("list")
                attrib = []
                for pt in self.points:
                    attrib.append(getattr(pt,name))
                return attrib
        else:
            print('attribute ' + name + ' does not exist')
            return False
    # 提取点群组 0（不在组内）和1（在组内）标记列表，用于c4d selection tag
    def GetPointGroupMark(self, name:str) -> list:
        groupmark = []
        if name in self.pointgroups:
            for pt in self.points:
                groupmark.append(getattr(pt,'group:' + name))
            return groupmark
        print('group ' + name + ' does not exist')
        return False
    # 提取点群组索引列表。未处理字符串包含于列表某一元素之中的情况
    def GetPointGroupIndex(self, groupname='') -> list:
        groupindex = []
        # 未传入参数，传入空字符串或传入空列表
        if len(groupname) == 0:
            return False
        # 传入字符串
        elif type(groupname) == str:
            if groupname in self.pointgroups:
                # print(groupname)
                for i in range(len(self.points)):
                    if getattr(self.points[i],'group:' + groupname) == 1:
                        groupindex.append(i)
                return groupindex
            else:
                print('group ' + groupname + ' does not exist')
                return False
        # 传入列表
        elif type(groupname) == list:
            for name in groupname:
                if name in self.pointgroups:
                    # print(name)
                    for i in range(len(self.points)):
                        if getattr(self.points[i],'group:' + name) == 1:
                            groupindex.append(i)
                else:
                    print('group '+ name + ' does not exist')
                    return False
            groupindex = list(set(groupindex))
            groupindex.sort()
            return groupindex
        else:
            print('invalid input name')
            return False        

    def GetPrimAttrib(self,name:str, index:int = -1):
        if name in self.primitiveattributes:
            if index in range(self.pointcount):
                attrib = getattr(self.primitives[index],name)
                return attrib
            else:
                attrib = []
                for prim in self.primitives:
                    attrib.append(getattr(prim,name))
                return attrib
        else:
            print('attribute ' + name + ' does not exist')
            return False
    
    def __repr__(self) -> str:
        output = "point number:\t"
        # point and pointgroup
        items = list(self.pointattributes) + self.pointgroups
        for item in items:
            output += item+"\t"
        output += "\n"
        for i in range(len(self.points)):
            output += str(i) + "\t"
            for item in self.pointattributes:
                ptattr = getattr(self.points[i],item)
                # print(item,str(print(getattr(pt,item))),getattr(pt,item))
                if self.pointattributes[item] != "str":
                    output += str(ptattr) + "\t"
                elif self.pointattributes[item] == "str":
                    if len(ptattr)==0:
                        output +=  "\t\t"
                    else:
                        output += ptattr + "\t"
            for item in self.pointgroups:
                gpattr = getattr(self.points[i],"group:"+item)
                output += str(gpattr) + "\t"
            output += "\n"
        # vertice
        output += "vertice number:\t"
        output += "point number\n"
        for i in range(len(self.primitives)):
            for j in range(len(self.primitives[i].vertices)):
                output += str(i) + ":" + str(j) + "\t"
                output += str(self.primitives[i].vertices[j].pointnumber) + "\n"
        output += "\n"
        # prim
        output += "primitive number:\t"
        items = list(self.primitiveattributes)
        for item in items:
            output += item+"\t"
        output += "\n"
        for i in range(len(self.primitives)):
            output += str(i) + "\t"
            for item in self.primitiveattributes:
                primattr = getattr(self.primitives[i],item)
                if self.primitiveattributes[item] != "str":
                    output += str(primattr) + "\t"
                elif self.primitiveattributes[item] == "str":
                    if len(primattr)==0:
                        output +=  "\t\t"
                    else:
                        output += primattr + "\t"
            output += "\n"
        return output

# All these functions modify the PDTGeo object---------------------------------------------------------------------------
class PDTFunction(object):
    @staticmethod
    def setpointgroup(geo:PDTGeo,name:str,point_num:int,value:int):
        if name in geo.pointgroups:
            setattr(geo.points[point_num],'group:' + name, value)
        else:
            geo.pointgroups.append(name)
            for i in range(len(geo.points)):
                if i == point_num:
                    setattr(geo.points[i],'group:' + name, value)
                else:
                    setattr(geo.points[i],'group:' + name, 0)

    @staticmethod
    def setpointattrib(geo:PDTGeo,name:str,pointnum:int,value):
        
        valuetype = type(value)
        if valuetype == np.ndarray:
            if str(value.shape) == "(3,)":
                if name not in geo.pointattributes:
                    geo.pointattributes[name] = "ndarray3"
                    for i in range(len(geo.points)):
                        if i == pointnum:
                            setattr(geo.points[i], name, value)
                        else:
                            setattr(geo.points[i], name, np.array([0.0,0.0,0.0]))
                else:
                    setattr(geo.points[pointnum], name, value)
        else:
            if name not in geo.pointattributes:
                geo.pointattributes[name] = valuetype
                for i in range(len(geo.points)):
                    if i == pointnum:
                        setattr(geo.points[i], name, value)
                    else:
                        setattr(geo.points[i], name, valuetype())
            else:
                setattr(geo.points[pointnum], name, value)
    @staticmethod
    def setprimattrib(geo:PDTGeo,name:str,primnum:int,value):
        valuetype = type(value)
        if valuetype == np.ndarray:
            if str(value.shape) == "(3,)":
                if name not in geo.primitiveattributes:
                    geo.primitiveattributes[name] = "ndarray3"
                    for i in range(len(geo.primitives)):
                        if i == primnum:
                            setattr(geo.primitives[i], name, value)
                        else:
                            setattr(geo.primitives[i], name, np.array([0.0,0.0,0.0]))
                else:
                    setattr(geo.primitives[primnum], name, value)
        else:
            if name not in geo.primitiveattributes:
                geo.primitiveattributes[name] = valuetype
                for i in range(len(geo.primitives)):
                    if i == primnum:
                        setattr(geo.primitives[i], name, value)
                    else:
                        setattr(geo.primitives[i], name, valuetype())
            else:
                setattr(geo.primitives[primnum], name, value)
        
    @staticmethod
    # 返回值是ptnum
    def addpoint(geo:PDTGeo,pos:np.ndarray)->int:
        point = PDTPoint(pos)
        for name, attrtype in geo.pointattributes.items():
            if name != "P":
                if attrtype == "ndarray3":
                    setattr(point,name,np.array([0.0,0.0,0.0]))
                else:
                    setattr(point,name,attrtype())
        for name in geo.pointgroups:
            setattr(point, 'group:' + name,0)
        geo.points.append(point)
        geo.pointcount += 1
        return geo.pointcount-1
    
    @staticmethod
    def addprim(geo:PDTGeo,primtype:str,ptnums:list)->int:
        geo.primitives.append(PDTPrimitive(ptnums,primtype))
        geo.primitivecount += 1
        return geo.primitivecount-1
    
    # 删除多个点按索引值遍历，输入顺序的ptnums，该函数从倒序处理删除;暂时未更新prim和vertice
    @staticmethod
    def removepoints(geo:PDTGeo,ptnums:list[int]):
        ptnums.reverse()
        for i in ptnums:
            geo.points.pop(i)
        geo.pointcount -= len(ptnums)

    @staticmethod
    def translate(geo:PDTGeo,t:np.ndarray):
        newgeo = copy.deepcopy(geo)
        # print(newgeo.primitiveattributes)
        for i in newgeo.points:
            i.P = i.P + t
        return newgeo
    
    @staticmethod
    def merge(geos:list[PDTGeo]):
        out = PDTGeo()
        for geo in geos:
            out.pointattributes.update(geo.pointattributes)
            out.primitiveattributes.update(geo.primitiveattributes)
            out.pointgroups.extend(geo.pointgroups)
            out.pointcount += geo.pointcount
        out.pointgroups = list(set(out.pointgroups))
        ptcount = 0
        for geo in geos:
            for pt in geo.points:
                for name,attrtype in out.pointattributes.items():
                    if name != "P":
                        if hasattr(pt,name) == False:
                            if attrtype == "ndarray3":
                                setattr(pt,name,np.array([0.0,0.0,0.0]))
                            else:
                                setattr(pt,name,attrtype())
                for name in out.pointgroups:
                    if hasattr(pt,'group:'+ name) == False:
                        setattr(pt,'group:'+ name,0)
                out.points.append(pt)
                ptcount += 1
            for prim in geo.primitives:
                for name,attrtype in out.primitiveattributes.items():
                    if hasattr(prim,name) == False:
                        if attrtype == "ndarray3":
                            setattr(prim,name,np.array([0.0,0.0,0.0]))
                        else:
                            setattr(prim,name,attrtype())
                # print("add prim")
                for vertice in prim.vertices:
                    vertice.pointnumber += ptcount - geo.pointcount
                out.primitives.append(prim)
        return out
# These functions return c4d data or geometry-----------------------------------------------------------------------
class C4DFunction(object):
    #convert float3 list to c4d.Vector list
    @staticmethod
    def ToC4dVector(f3list:list) -> list[c4d.Vector]:
        c4dpos:list[c4d.Vector] = []
        for i in f3list:
            c4dpos.append(c4d.Vector(i[0],i[1],i[2]))
        return c4dpos

    #generate point cloud mesh from c4d.Vector pos list
    @staticmethod
    def PostoMesh(pos) -> c4d.BaseObject:
        poly = c4d.BaseObject(5100)
        poly.ResizeObject(len(pos))
        poly.SetAllPoints(pos)
        poly.Message(c4d.MSG_UPDATE)
        return poly
    # get a position list from a polygon object based on selection tag name
    @staticmethod
    def GroupPos(geo:c4d.PolygonObject = c4d.BaseObject(5100), group:str = '') -> list[c4d.Vector]:
        pts = geo.GetAllPoints()
        tags = geo.GetTags()
        pts_group = []
        for i in range(len(tags)):
            if type(tags[i]) == c4d.SelectionTag and tags[i][c4d.ID_BASELIST_NAME] == group:
                sel = tags[i].GetBaseSelect()
                marklist = sel.GetAll(len(pts))
                for j in range(len(marklist)):
                    if marklist[j] == 1:
                        pts_group.append(pts[j])
                return pts_group
        print('selection tag '+ group + ' does not exist')
        return False
    # get position list from point index list
    @staticmethod
    def IndexPos(geo:c4d.PolygonObject, index:list) -> list[c4d.Vector]:
        pos:list[c4d.Vector] = []
        pts = geo.GetAllPoints()
        for i in index: 
            pos.append(pts[i])
        return pos
    
    @staticmethod
    def AddCloner(linkobj:c4d.BaseObject, name:str = "Cloner", groupname:str = [], color = c4d.Vector(1)) -> c4d.BaseObject:
        cln = c4d.BaseObject(1018544)
        cln[c4d.ID_BASELIST_NAME] = name
        cln[c4d.MGCLONER_VOLUMEINSTANCES_MODE] = 1
        cln[c4d.ID_MG_MOTIONGENERATOR_MODE] = 0
        cln[c4d.MG_OBJECT_LINK] = linkobj
        cln[c4d.MG_POLY_MODE_] = 0
        cln[c4d.MG_POLY_SELECTION] = groupname
        cln[c4d.ID_MG_TRANSFORM_COLOR] = color
        return cln
    
    @staticmethod
    def AddSphere(name:str = "Sphere", rad:float = 1, seg:int = 20) ->c4d.BaseObject:
        sph = c4d.BaseObject(5160)
        sph[c4d.ID_BASELIST_NAME] = name
        sph[c4d.PRIM_SPHERE_RAD] = rad
        sph[c4d.PRIM_SPHERE_TYPE] = 4
        sph[c4d.PRIM_SPHERE_SUB] = 10
        sphphone = c4d.BaseTag(5612)
        sphphone[c4d.PHONGTAG_PHONG_ANGLELIMIT] = 1
        sph.InsertTag(sphphone)
        return sph
    
    @staticmethod
    def CreateMesh(geo:PDTGeo,targetpoint:np.ndarray = np.empty(shape=(0,3))) ->c4d.BaseObject:
        mesh = c4d.BaseObject(5100)
        if targetpoint.size == 0:
            pos = geo.GetP()
        elif targetpoint.ndim == 1:
            pos = []
            for i in targetpoint:
                pos.append(PDTTranslation.nptoc4d(geo.points[i].P))
        else:
            print("targetpoint error, use 1 dim np.ndarray")
            return None

        mesh.ResizeObject(len(pos))
        mesh.SetAllPoints(pos)
        mesh.Message(c4d.MSG_UPDATE)

        if len(geo.pointgroups) == 0:
            return mesh
        else:
            for name in geo.pointgroups:
                #creat point selection tags
                tsel = c4d.BaseTag(c4d.Tpointselection)
                tsel[c4d.ID_BASELIST_NAME] = name
                #the returned BaseSelect obj can be modified and automatically updated back to selection tag
                sel = tsel.GetBaseSelect()
                marklist = []
                if targetpoint.size == 0:
                    for i in range(len(pos)):
                        if getattr(geo.points[i],'group:'+name) == 1:
                            marklist.append(1)
                        else:
                            marklist.append(0)
                else:
                    for i in targetpoint:
                        if getattr(geo.points[i],'group:'+name) == 1:
                            marklist.append(1)
                        else:
                            marklist.append(0)
                # print(marklist)
                #select points from mark list, 1 is selected, 0 is unselected
                sel.SetAll(marklist)
                #attach tag to mesh obj
                mesh.InsertTag(tsel)
            return mesh

        
    @staticmethod
    def CreateAtom(geo:PDTGeo, attrname:str="label", targetpoint:np.ndarray = np.empty(shape=(0,3))):
        rootatom = c4d.BaseObject(5140)
        rootatom[c4d.ID_BASELIST_NAME] = 'atom'
        rootclone = c4d.BaseObject(5140)
        rootclone[c4d.ID_BASELIST_NAME] = 'clone'
        rootclone.InsertUnder(rootatom)

        mesh = C4DFunction.CreateMesh(geo,targetpoint)
        mesh.InsertUnder(rootatom)
        labels = geo.GetPointAttrib(attrname)
        labels = list(set(labels))
        for label in labels:
            
            index = geo.GetPointGroupIndex(label)[0]
            cd = PDTTranslation.nptoc4d(geo.GetPointAttrib("Cd",index))
            pscale = geo.GetPointAttrib("pscale",index)

            sph = C4DFunction.AddSphere(label,pscale)
            cln = C4DFunction.AddCloner(mesh, label, label, cd)

            sph.InsertUnder(cln)
            cln.InsertUnder(rootclone)
        return rootatom

    @staticmethod
    def CreateBond(geo_bond:PDTGeo, attrname:str="label", targetprim:np.ndarray = np.empty(shape=(0,3))) ->c4d.BaseObject:
        if attrname in geo_bond.primitiveattributes.keys():
            rootsweep = c4d.BaseObject(5140)
            names = set(geo_bond.GetPrimAttrib(attrname))
            # print(names)
            for name in names:
                spl = c4d.BaseObject(5101)
                linepos = []
                cd = c4d.Vector(0)
                if targetprim.size == 0:
                    for prim in geo_bond.primitives:
                        if getattr(prim,attrname) == name:
                            cd = getattr(prim,"Cd")
                            for vertex in prim.vertices:
                                linepos.append(geo_bond.points[vertex.pointnumber].P)
                else:
                    for i in targetprim:
                        if getattr(geo_bond.primitives[i],attrname) == name:
                            cd = getattr(geo_bond.primitives[i],"Cd")
                            for vertex in geo_bond.primitives[i].vertices:
                                linepos.append(geo_bond.points[vertex.pointnumber].P)
                #set point count and segment count
                spl.ResizeObject(len(linepos), int(len(linepos)/2))

                #set point positions and assign point to segment
                # print(linepos)

                spl.SetAllPoints(PDTTranslation.listtoc4d(linepos))
                for j in range(int(len(linepos)/2)):
                    spl.SetSegment(j, 2, False)
                
                #setup sweep
                spl[c4d.ID_BASELIST_NAME] = name
                spl.Message(c4d.MSG_UPDATE)
                circle = c4d.BaseObject(5181)
                circle[c4d.PRIM_CIRCLE_RADIUS] = 10
                sweep = c4d.BaseObject(5118)
                sweepphone = c4d.BaseTag(5612)
                sweepphone[c4d.PHONGTAG_PHONG_ANGLELIMIT] = 1
                sweep.InsertTag(sweepphone)
                sweep[c4d.ID_BASEOBJECT_USECOLOR] = 2
                sweep[c4d.ID_BASEOBJECT_COLOR] = PDTTranslation.nptoc4d(cd)
                sweep[c4d.ID_BASELIST_NAME] = name
                spl.InsertUnder(sweep)
                circle.InsertUnder(sweep)
                sweep.InsertUnder(rootsweep)
            return rootsweep
        else:
            print("PDTGeo doesnot has '" + attrname + "' prim attribute")
            return False
        
    @staticmethod
    def Vis_pt(pts:list[c4d.Vector]):
        doc:c4d.documents.BaseDocument
        ins = c4d.BaseObject(5126)
        ins[c4d.INSTANCEOBJECT_RENDERINSTANCE_MODE] = 2
        ins[c4d.INSTANCEOBJECT_DRAW_MODE] = 3
        mtx = []
        for pt in pts:
            mtx.append(c4d.Matrix(off=pt))
        ins.SetInstanceMatrices(mtx)
        doc.InsertObject(ins)
        return True
    
def main():
    pts = np.arange(12).reshape(4,3)
    # print(pts)
    # for i in pts:
    #     print(i)
    a = PDTGeo(pts)
    PDTFunction.setpointgroup(a,'group1', 0, 1)
    PDTFunction.setpointattrib(a,'Cd',1,np.array((0.1,0.1,0.1)))
    PDTFunction.setpointattrib(a,'pscale',1,1.0)
    PDTFunction.setpointattrib(a,'path',2,"obj")
    PDTFunction.addpoint(a,np.array([0.4,0.0,0.0]))
    PDTFunction.addprim(a,"polyline",[0,1])
    PDTFunction.setprimattrib(a,"Cd",0,np.array((1.0,0.0,0.0)))
    # print(a.pointattributes)
    # print(a.pointattributes)
    # b = PDTGeo([c4d.Vector()]*4)
    # PDTFunction.setpointgroup(b,'group2', 0, 1)
    # PDTFunction.addprim(b,'polyline',[1,2])
    # PDTFunction.addprim(b,'polyline',[3,4])
    # PDTFunction.setprimattrib(b,'symbol',1,'C')
    # # print(a)
    
    c = PDTFunction.translate(a,np.array((1,0,0)))
    geo = []
    for i in range(2):
        geo.append(PDTFunction.translate(c,np.array((1.0,0.0,0.0))))
    # print(geo)
    d = PDTFunction.merge(geo)
    print(a,c,d)
    # print(d.primitives[0].pointcount)
    # print(b)
    # geo = PDTFunction.merge([a,b])
    # print(geo)
    # PDTFunction.setprimattrib(a,'symbol',0,'D')
    # geo = PDTFunction.merge([a,b])
    # for i in range(a.pointcount):
    #     PDTFunction.addpoint(a,c4d.Vector(0))
    #     PDTFunction.setpointattrib(a,"pscale",i,0.1)
    #     print(a.points[i].P,getattr(a.points[i],'group:group1'),getattr(a.points[i],'group:group2'),a.points[i].Cd,a.points[i].pscale,a.points[i].path)
    # print(geo.primitivecount)
    # for i in range(a.primitivecount):
    #     print(a.primitives[i].symbol)
    # print(a.pointattributes)
    # print(a.pointcount)
if __name__ == '__main__':
    main()