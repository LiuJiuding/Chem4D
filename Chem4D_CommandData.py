import c4d
import typing
import os
import gemmi
import numpy as np
from PDT import PDTGeo, PDTFunction, C4DFunction
from rdkit import Chem

from c4d import BaseContainer
PLUGIN_ID = 1059606
PLUGI_NNAME = "Chem4DCommand"
VERSION = "1.1.0"

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
        # return mol

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
    
    @staticmethod
    def ReadCif(path, a_min=0, a_max=1, b_min=0, b_max=1, c_min=0, c_max=1):
        #read cif file
        cif_doc = gemmi.cif.read(path)
        LCO = gemmi.make_small_structure_from_block(cif_doc.sole_block())
        # sg = LCO.spacegroup_hm
        cell = LCO.cell
        op_list = cif_doc[0].find_values('_symmetry_equiv_pos_as_xyz')
        replace_dict = {"'":"",".3333":"1/3",".6666":"2/3",".5":"1/2",".75":"3/4",".25":"1/4"}
        constmat = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]]

        geo = PDTGeo()

        # fullpos = []
        #apply symmetry operations to element sites
        
        for site in LCO.sites:
            # print(site)
            poslist = []
            #get operator list and apply op
            for i in range(len(op_list)):
                newop = op_list[i]
                
                #replace unreadable strings
                for j in replace_dict:
                    newop = newop.replace(j,replace_dict[j])
                #get gemmi symmetry operator
                op = gemmi.Op(newop)
                fracpos = op.apply_to_xyz([site.fract[0],site.fract[1],site.fract[2]])
                for trans in constmat:
                    transpos = [fracpos[0]+trans[0],fracpos[1]+trans[1],fracpos[2]+trans[2]]
                    poslist.append(transpos)

            #delete repeated position
            poslist = np.unique(poslist, axis=0).tolist()

            #delete position outside unit cell
            min = -0.01
            max = 1.01
            unitpos = []
            for i in poslist:
                if i[0]<max and i[1]<max and i[2]<max and i[0]>min and i[1]>min and i[2]>min:
                    unitpos.append(i)
            poslist.clear()

            #expand unit cell
            pos = []
            for i in range(a_min,a_max,1):
                for j in range(b_min,b_max,1):
                    for k in range(c_min,c_max,1):
                        for x in unitpos:
                            pos.append([x[0]+i,x[1]+j,x[2]+k])   

            #delete repeated position
            pos = np.unique(pos, axis=0).tolist()

            #transform to orth pos
            for i in range(len(pos)):
                pos[i] = c4d.Vector(cell.orthogonalize(gemmi.Fractional(pos[i][0],pos[i][1],pos[i][2]))[0],cell.orthogonalize(gemmi.Fractional(pos[i][0],pos[i][1],pos[i][2]))[1],cell.orthogonalize(gemmi.Fractional(pos[i][0],pos[i][1],pos[i][2]))[2])

            #creat atom at position
            for i in pos:
                ptnum = PDTFunction.addpoint(geo, i*10)
                PDTFunction.setpointattrib(geo, "symbol", ptnum, site.type_symbol)
                PDTFunction.setpointgroup(geo, site.type_symbol, ptnum, 1 )

        return geo
class ChemDialog(c4d.gui.GeDialog, CHEM4DHelper):
    ID_GROUP = 100
    ID_FILEPATH = 1000
    ID_BUTTON_CREATE = 1001
    ID_OBJECT_LINK = 1002
    ID_UPDATE = 1003
    # The IDs of link boxes and their accepted node types as a hashmap.
    LNK_ITEMS: dict[int, int] = {
        1002: c4d.Onull # This link box accepts any node of type Onull
    }
    def __init__(self) -> None:
        self._linkCache: dict[int, c4d.BaseList2D] = {}
        super().__init__()

    def InitValues(self):
        #read property file
        directory, _ = os.path.split(__file__)
        fn = os.path.join(directory, "AtomProperties.txt")
        self.read_dprop = self.ReadProperty(fn)

    def _validateLink(self, gid: int) -> bool:
        """Validates the link box gadget at #gid.

        When #gid has an illegal value, its cached state will be restored.
        """
        # The linked item is valid, we have nothing to do. 
        node: typing.Optional[c4d.BaseList2D] = self._getLink(gid)
        if node is not None:
            return True
        
        # It is not, but the cache somehow went out of whack.
        if gid not in self._linkCache:
            raise RuntimeError(f"Validation event for uncached parameter value for gid: {gid}")
        
        # Restore the previous state.
        res: bool = self._setLink(gid, self._linkCache[gid])
        if res:
            del(self._linkCache[gid])
        
        return res
        
    # filter target object
    def _getLink(self, gid: int) -> typing.Optional[c4d.BaseList2D]:
        """Gets the node value of the link box at #gid.

        This will use LinkDialog.LNK_ITEMS to determine the type of nodes #gid is supposed to be
        able to hold. When #gid holds an 'illegal' type, this will return None.
        """
        if gid not in ChemDialog.LNK_ITEMS.keys():
            return None

        gadget: c4d.gui.LinkBoxGui = self.FindCustomGui(gid, c4d.CUSTOMGUI_LINKBOX)
        if not isinstance(gadget, c4d.gui.LinkBoxGui):
            raise RuntimeError(f"{gadget = }, {gid = }")
        
        # Get the accepted node types from LNK_ITEMS and rely on the type checking of .GetLink. We
        # could also just retrieve the node without that argument and check ourselves when we want
        # to do more fancy things like (Obase | Mbase), i.e., accept a BaseObject or a BaseMaterial.
        typeId: int = ChemDialog.LNK_ITEMS[gid]
        nullobj = gadget.GetLink(c4d.documents.GetActiveDocument(), typeId)
        if nullobj == None:
            return None
        id_bc_list = nullobj.GetUserDataContainer()
        try:
            name = id_bc_list[0][1].GetString(1)
        except IndexError:
            return None
        if name == "molecule":
            return nullobj
    def _setLink(self, gid: int, value: c4d.BaseList2D) -> bool:
        """Sets the node value of the link box at #gid.
        """
        if gid not in ChemDialog.LNK_ITEMS.keys():
            return False
        
        gadget: c4d.gui.LinkBoxGui = self.FindCustomGui(gid, c4d.CUSTOMGUI_LINKBOX)
        if not isinstance(gadget, c4d.gui.LinkBoxGui):
            raise RuntimeError(f"{gadget = }, {gid = }")
        
        return gadget.SetLink(value)
    def CreateLayout(self):
        """This Method is called automatically when Cinema 4D Create the Layout (display) of the Dialog."""
        # Defines the title of the Dialog
        self.SetTitle(PLUGI_NNAME+VERSION)


        self.GroupBegin(self.ID_GROUP, c4d.BFH_SCALEFIT | c4d.BFV_SCALEFIT, cols= 2, rows= 2, title= "", groupflags= 0, initw= 0, inith= 0)
        
        # set gui by following line
        # customdata[c4d.FILENAME_SAVE] = True
        self.AddCustomGui(self.ID_FILEPATH, c4d.CUSTOMGUI_FILENAME, name= "File", flags= c4d.BFH_SCALEFIT, minw= 10, minh= 10, customdata= c4d.BaseContainer())
        self.AddButton(self.ID_BUTTON_CREATE, c4d.BFH_FIT, initw= 50, inith= 10, name= "Create")
        self.AddCustomGui(self.ID_OBJECT_LINK, c4d.CUSTOMGUI_LINKBOX, name= "Object", flags= c4d.BFH_SCALEFIT, minw= 10, minh= 10, customdata= c4d.BaseContainer())
        self.AddButton(self.ID_UPDATE, c4d.BFH_FIT, initw= 50, inith= 10, name= "Update")
        self.GroupEnd()
        return True
    
    def Command(self, mid: int, msg: BaseContainer) -> bool:
        
        """This Method is called automatically when the user clicks on a gadget and/or changes its value this function will be called.
        It is also called when a string menu item is selected.

        Args:
            messageId (int): The ID of the gadget that triggered the event.
            bc (c4d.BaseContainer): The original message container.

        Returns:
            bool: False if there was an error, otherwise True.
        """
        # One of the link boxes has been set, we validate the new state and optionally revert to
        # the old state. This could also be done with BFM_ACTION in Message().
        if mid == self.ID_OBJECT_LINK:
            self._validateLink(mid)

        # User click on Creat
        if mid == self.ID_BUTTON_CREATE:
            #check path and creat rdkit mol
            path = self.GetFilename(self.ID_FILEPATH)
            if path.endswith(".mol") or path.endswith(".mol2"):
                print(self.GetFilename(self.ID_FILEPATH))
                mol = self.ReadMol(self.GetFilename(self.ID_FILEPATH))
                if mol == False:
                    return False
                # print(self.GetFilename(self.ID_FILEPATH), "Create MOL")
                geo = mol[0]
                geo_bond = mol[1]
                parameters = []
                for i in geo.pointgroups:
                    parameters.append(self.read_dprop[i][0])
                    parameters.append(self.read_dprop[i][1])

                doc = c4d.documents.GetActiveDocument()

                rootobj = c4d.BaseObject(5140)
                rootobj[c4d.ID_BASELIST_NAME] = os.path.basename(self.GetFilename(self.ID_FILEPATH))
                rootatom = c4d.BaseObject(5140)
                rootatom[c4d.ID_BASELIST_NAME] = 'atom'
                rootbond = c4d.BaseObject(5140)
                rootbond[c4d.ID_BASELIST_NAME] = 'bond'
                rootatom.InsertUnder(rootobj)
                rootbond.InsertUnder(rootobj)
                
                pmesh = C4DFunction.CreateMesh(geo)
                # pos = geo.GetP()
                # pmesh = C4DFunction.PostoMesh(pos)
                pmesh[c4d.ID_BASELIST_NAME] = "PositionMesh"
                pmesh.InsertUnder(rootobj)

                # create control parameters (user data) for root obj------------------------------------
                # structure type
                bc_type = c4d.GetCustomDataTypeDefault(c4d.DTYPE_STRING)
                bc_type[c4d.DESC_NAME] = "molecule"
                bc_type[c4d.DESC_HIDE] = True
                rootobj.AddUserData(bc_type)
                # group: general parameters
                bc_group_general = c4d.GetCustomDataTypeDefault(c4d.DTYPE_GROUP)
                bc_group_general[c4d.DESC_NAME] = "General"
                bc_group_general[c4d.DESC_TITLEBAR ] = True
                descId_group_general = rootobj.AddUserData(bc_group_general)
                # group: atom parameters
                bc_group_atom = c4d.GetCustomDataTypeDefault(c4d.DTYPE_GROUP)
                bc_group_atom[c4d.DESC_NAME] = "Atom"
                bc_group_atom[c4d.DESC_TITLEBAR ] = True
                descId_group_atom = rootobj.AddUserData(bc_group_atom)
                # build type, ball ,stick or both
                bc_buildtype = c4d.GetCustomDataTypeDefault(c4d.DTYPE_LONG)
                bc_buildtype[c4d.DESC_NAME] = "Build type"
                bc_buildtype[c4d.DESC_PARENTGROUP ] = descId_group_general
                bc_buildtype[c4d.DESC_CUSTOMGUI] = 200000281
                bc = c4d.BaseContainer()
                bc[0] = "Ball"
                bc[1] = "Stick"
                bc[2] = "B&S"
                bc_buildtype[c4d.DESC_CYCLE] = bc
                rootobj.AddUserData(bc_buildtype)
                # uniform overall scale of atoms
                bc_uniscale = c4d.GetCustomDataTypeDefault(c4d.DTYPE_REAL)
                bc_uniscale[c4d.DESC_NAME] = "Uniform Scale"
                bc_uniscale[c4d.DESC_PARENTGROUP ] = descId_group_general
                descId_uniscale = rootobj.AddUserData(bc_uniscale)
                rootobj[descId_uniscale] = 1
                # segment of sphere
                bc_segment = c4d.GetCustomDataTypeDefault(c4d.DTYPE_LONG)
                bc_segment[c4d.DESC_NAME] = "Sphere Segment"
                bc_segment[c4d.DESC_PARENTGROUP] = descId_group_general
                bc_segment[c4d.DESC_MIN] = 0
                descId_segment = rootobj.AddUserData(bc_segment)
                rootobj[descId_segment] = 10            

                #iterate atom symbols and build hierarchy--------------------------------------------------
                indnum = 0
                rad = 1

                for name in geo.pointgroups:
                    # create userdata on rootobj-------------------------------------
                    # use template to set type and name
                    bc_radius = c4d.GetCustomDataTypeDefault(c4d.DTYPE_REAL)
                    bc_radius[c4d.DESC_NAME] = "Radius " + name
                    bc_radius[c4d.DESC_MIN] = 0
                    bc_radius[c4d.DESC_MAX] = 200
                    bc_radius[c4d.DESC_PARENTGROUP ] = descId_group_atom
                    # get id
                    descId_radius = rootobj.AddUserData(bc_radius)
                    # set value
                    rootobj[descId_radius] = parameters[2*indnum]*rad

                    bc_color = c4d.GetCustomDataTypeDefault(c4d.DTYPE_COLOR)
                    bc_color[c4d.DESC_NAME] = "Color " + name
                    bc_color[c4d.DESC_PARENTGROUP ] = descId_group_atom
                    descId_color = rootobj.AddUserData(bc_color)
                    rootobj[descId_color] = parameters[2*indnum+1]

                    #creat atoms------------------------------------------------------
                    #setup cloner
                    sph = C4DFunction.AddSphere(name)
                    cln = C4DFunction.AddCloner(pmesh,name)

                    cln[c4d.ID_MG_TRANSFORM_SCALE] = c4d.Vector(rootobj[descId_radius])
                    cln[c4d.ID_MG_TRANSFORM_COLOR] = rootobj[descId_color]

                    sph.InsertUnder(cln)
                    cln.InsertUnderLast(rootatom)

                    #creat bond splines-----------------------------------------------
                    sweep = C4DFunction.CreateSweep(geo_bond,name)
                    sweep[c4d.ID_BASEOBJECT_COLOR] = rootobj[descId_color]
                    sweep.InsertUnderLast(rootbond)

                    indnum = indnum +1

                doc.InsertObject(rootobj)
                c4d.EventAdd()
            if path.endswith(".cif"):
                pass
        return super().Command(mid, msg)
    def Message(self, msg: BaseContainer, result: BaseContainer) -> int:
        # This is a drag and drop event for a link box, we store its current state before the drag
        # and drop has been carried out.
        if msg.GetId() == c4d.MSG_DESCRIPTION_CHECKDRAGANDDROP:
            gid: int = msg.GetLong(c4d.LINKBOX_ACCEPT_MESSAGE_CONTROL_ID)
            if gid == self.ID_OBJECT_LINK:
                # print("id in target")
                self._linkCache[gid] = self._getLink(gid)
                # print(self._getLink(gid))
        return super().Message(msg, result)
    
class ChemCommand(c4d.plugins.CommandData):
    """Command Data class that holds the ExampleDialog instance."""
    dialog = None
    
    def Execute(self, doc):
        """Called when the user executes a command via either CallCommand() or a click on the Command from the extension menu.

        Args:
            doc (c4d.documents.BaseDocument): The current active document.

        Returns:
            bool: True if the command success.
        """
        # Creates the dialog if its not already exists
        if self.dialog is None:
            self.dialog = ChemDialog()

        # Opens the dialog
        return self.dialog.Open(dlgtype=c4d.DLG_TYPE_ASYNC, pluginid=PLUGIN_ID, defaultw=400, defaulth=32)

    def RestoreLayout(self, sec_ref):
        """Used to restore an asynchronous dialog that has been placed in the users layout.

        Args:
            sec_ref (PyCObject): The data that needs to be passed to the dialog.

        Returns:
            bool: True if the restore success
        """
        # Creates the dialog if its not already exists
        if self.dialog is None:
            self.dialog = ChemDialog()

        # Restores the layout
        return self.dialog.Restore(pluginid=PLUGIN_ID, secret=sec_ref)
    
# main
if __name__ == "__main__":
    # Registers the plugin
    c4d.plugins.RegisterCommandPlugin(id=PLUGIN_ID,
                                      str="Chem4D",
                                      info=0,
                                      help="Display a basic GUI",
                                      dat=ChemCommand(),
                                      icon=None)

    

#             t[c4d.TPYTHON_CODE] = '''
# from typing import Optional
# import c4d

# doc: c4d.documents.BaseDocument # The document evaluating this tag
# op: c4d.BaseTag # The Python scripting tag
# flags: int # c4d.EXECUTIONFLAGS
# priority: int # c4d.EXECUTIONPRIORITY
# tp: Optional[c4d.modules.thinkingparticles.TP_MasterSystem] # Particle system
# bt: Optional[c4d.threading.BaseThread] # The thread executing this tag

# def main() -> None:
#     # Called when the tag is executed. It can be called multiple time per frame. Similar to TagData.Execute.
#     # Write your code here
#     root = op.GetObject()
#     root_bond = root.GetDownLast().GetPred().GetChildren()
#     root_atom = root.GetDownLast().GetChildren()
#     # print(root_bond,root_atom)
#     t = root.GetTag(c4d.Tpython)
#     for i in range(len(root_atom)):
#         root_bond[i][c4d.ID_BASEOBJECT_COLOR] = t[c4d.ID_USERDATA, 2*i+2]
#         root_atom[i][c4d.ID_MG_TRANSFORM_SCALE] = c4d.Vector(t[c4d.ID_USERDATA, 2*i+1])
#         root_atom[i][c4d.ID_MG_TRANSFORM_COLOR] = t[c4d.ID_USERDATA, 2*i+2]

# '''