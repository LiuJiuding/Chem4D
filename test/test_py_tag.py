import c4d
import copy
from typing import Optional
import numpy as np
doc: c4d.documents.BaseDocument # The document evaluating this tag
op: c4d.BaseTag # The Python scripting tag
flags: int # c4d.EXECUTIONFLAGS
priority: int # c4d.EXECUTIONPRIORITY
tp: Optional[c4d.modules.thinkingparticles.TP_MasterSystem] # Particle system
bt: Optional[c4d.threading.BaseThread] # The thread executing this tag


# root = doc.GetActiveObject()
# root_bond = root.GetDownLast().GetPred().GetChildren()
# root_atom = root.GetDownLast().GetChildren()
# print(root_bond,root_atom)
# t = root.GetTag(c4d.Tpython)
# for i in range(len(root_atom)):
#     root_bond[i][c4d.ID_BASEOBJECT_COLOR] = t[c4d.ID_USERDATA, 2*i+2]
#     root_atom[i][c4d.ID_MG_TRANSFORM_SCALE] = c4d.Vector(t[c4d.ID_USERDATA, 2*i+1])
#     root_atom[i][c4d.ID_MG_TRANSFORM_COLOR] = t[c4d.ID_USERDATA, 2*i+2]
"""
def message(id: int, data: object) -> bool:
    # Called when the tag receives messages. Similar to TagData.Message.
    # Write your code here
    return super().Message(id, data)

def draw(bd: c4d.BaseDraw) -> bool:
    # Called to display some visual element in the viewport. Similar to TagData.Draw.
    # Write your code here
    return True
"""
class Test(object):
    def __init__(self,name,v) -> None:
        self.name = name
        self.v = v
    def convert(self):
        return c4d.Vector(float(self.v[0]),float(self.v[1]),float(self.v[2]))
def main():
    a = Test("a",np.array([1,1,1]))
    b = copy.deepcopy(a)
    a.v = np.array([1,2,1])
    # print(a.v,a.v[0],b.v)
    print(a.convert(),b.convert())

if __name__ == '__main__':
    main()