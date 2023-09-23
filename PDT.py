############ PDT v1.1.0 ###############
import c4d

class PDTPoint(object):
    def __init__(self,pos:c4d.Vector = c4d.Vector(0,0,0) ) -> None:
        self.P = pos
class PDTVertice(object):
    def __init__(self, ptnum:int) -> None:
        self.pointnumber = ptnum
class PDTPrimitive(object):
    def __init__(self,ptnums:list) -> None:
        self.primtype = "polyline"
        self.pointcount = len(ptnums)
        self.vertices = []
        for i in ptnums:
            v = PDTVertice(i)
            self.vertices.append(v)
class PDTGeo(object):
    def __init__(self, pts:list = []) -> None:
        self.points = []
        for i in range(len(pts)):
            self.points.append(PDTPoint(pts[i]))        
        self.primitives = []
        self.vertices = []

        self.pointcount = len(pts)
        self.primitivecount = len(self.primitives)

        self.pointattributes = {"P":type(c4d.Vector())}
        self.primitiveattributes = {}

        self.pointgroups = []

    def GetP(self,names:list = []) -> list:
        pos_list = []
        if len(names) == 0:
            for pt in self.points:
                pos_list.append(pt.P)
            return pos_list
        for i in range(len(names)):
            if names[i] in self.pointgroups:
                for pt in self.points:
                    if getattr(pt,'group:'+names[i]) == 1:
                        pos_list.append(pt.P)
            else:
                print('group '+ names[i] + ' does not exist')
                return False
        return pos_list
    def GetPointAttrib(self, name = str):
        attrib = []
        if name in self.pointattributes:
            for pt in self.points:
                attrib.append(getattr(pt,name))
            return attrib
        print('attribute ' + name + ' does not exist')
        return False
    def GetPointGroup(self, name:str) -> list:
        groupattrib = []
        if name in self.pointgroups:
            for pt in self.points:
                groupattrib.append(getattr(pt,'group:' + name))
            return groupattrib
        print('group ' + name + ' does not exist')
        return False
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
    def addpoint(geo:PDTGeo,pos:c4d.Vector)->int:
        point = PDTPoint(pos)
        for name in geo.pointattributes:
            if name != "P":
                setattr(point,name,geo.pointattributes[name]())
        for name in geo.pointgroups:
            setattr(point, 'group:' + name,0)
        geo.points.append(point)
        geo.pointcount += 1
        return geo.pointcount-1
    @staticmethod
    def addprim(geo:PDTGeo,primtype:str,ptnums:list)->int:
        if primtype == "polyline":
            geo.primitives.append(PDTPrimitive(ptnums))
            geo.primitivecount += 1
        return geo.primitivecount-1
class C4DFunction(object):
    #convert float3 list to c4d.Vector list
    @staticmethod
    def toc4dvector(f3list:list):
        c4dpos = []
        for i in f3list:
            c4dpos.append(c4d.Vector(i[0],i[1],i[2]))
        return c4dpos

    #generate point cloud mesh from c4d.Vector pos list
    @staticmethod
    def postomesh(pos):
        poly = c4d.BaseObject(5100)
        poly.ResizeObject(len(pos))
        poly.SetAllPoints(pos)
        poly.Message(c4d.MSG_UPDATE)
        return poly

def main():
    a = PDTGeo([c4d.Vector()]*4)
    PDTFunction.setpointgroup(a,'group1', 1, 1)
    PDTFunction.setpointattrib(a,'Cd',0,c4d.Vector(0.1,0.1,0.1))
    PDTFunction.setpointattrib(a,'path',2,"obj")
    PDTFunction.addpoint(a,c4d.Vector(0))
    PDTFunction.setpointgroup(a,'group2', 4, 1)
    PDTFunction.addprim(a,'polyline',[1,2])
    PDTFunction.addprim(a,'polyline',[3,4])
    PDTFunction.setprimattrib(a,'symbol',1,'C')
    PDTFunction.setprimattrib(a,'symbol',0,'D')
    # for i in range(a.pointcount):
    #     PDTFunction.addpoint(a,c4d.Vector(0))
    #     PDTFunction.setpointattrib(a,"pscale",i,0.1)
    #     print(a.points[i].P,getattr(a.points[i],'group:group1'),getattr(a.points[i],'group:group2'),a.points[i].Cd,a.points[i].pscale,a.points[i].path)
    print(a.primitivecount)
    for i in range(a.primitivecount):
        print(a.primitives[i].symbol)
    # print(a.pointattributes)
    # print(a.pointcount)
if __name__ == '__main__':
    main()