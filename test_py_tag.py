import c4d
from typing import Optional

doc: c4d.documents.BaseDocument # The document evaluating this tag
op: c4d.BaseTag # The Python scripting tag
flags: int # c4d.EXECUTIONFLAGS
priority: int # c4d.EXECUTIONPRIORITY
tp: Optional[c4d.modules.thinkingparticles.TP_MasterSystem] # Particle system
bt: Optional[c4d.threading.BaseThread] # The thread executing this tag


root = doc.GetActiveObject()
root_bond = root.GetDownLast().GetPred().GetChildren()
root_atom = root.GetDownLast().GetChildren()
print(root_bond,root_atom)
t = root.GetTag(c4d.Tpython)
for i in range(len(root_atom)):
    root_bond[i][c4d.ID_BASEOBJECT_COLOR] = t[c4d.ID_USERDATA, 2*i+2]
    root_atom[i][c4d.ID_MG_TRANSFORM_SCALE] = c4d.Vector(t[c4d.ID_USERDATA, 2*i+1])
    root_atom[i][c4d.ID_MG_TRANSFORM_COLOR] = t[c4d.ID_USERDATA, 2*i+2]
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