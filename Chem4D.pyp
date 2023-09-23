import c4d
import os
from PDT import PDTGeo, PDTFunction, C4DFunction
from rdkit import Chem
PLUGIN_ID = 1059606
PLUGI_NNAME = "Chem4D"
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
        return mol
    @staticmethod
    def BuildPDTGeo(mol):
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

class ChemDialog(c4d.gui.GeDialog, CHEM4DHelper):
    ID_GROUP = 100
    ID_FILEPATH = 1000
    ID_BUTTON_CREATE = 1001

    def InitValues(self):
        #read property file
        directory, _ = os.path.split(__file__)
        fn = os.path.join(directory, "AtomProperties.txt")
        self.read_dprop = self.ReadProperty(fn)

    def CreateLayout(self):
        """This Method is called automatically when Cinema 4D Create the Layout (display) of the Dialog."""
        # Defines the title of the Dialog
        self.SetTitle(PLUGI_NNAME+VERSION)

        # Creates a Ok and Cancel Button
        self.GroupBegin(self.ID_GROUP, c4d.BFH_SCALEFIT | c4d.BFV_SCALEFIT, cols= 2, rows= 1, title= "", groupflags= 0, initw= 0, inith= 0)
        customdata = c4d.BaseContainer()
        # set gui by following line
        # customdata[c4d.FILENAME_SAVE] = True
        self.AddCustomGui(self.ID_FILEPATH, c4d.CUSTOMGUI_FILENAME, name= "File", flags= c4d.BFH_SCALEFIT, minw= 10, minh= 10, customdata= customdata)
        self.AddButton(self.ID_BUTTON_CREATE, c4d.BFH_FIT, initw= 50, inith= 10, name= "Create")
        self.GroupEnd()
        return True
    
    def Command(self, messageId, bc):
        """This Method is called automatically when the user clicks on a gadget and/or changes its value this function will be called.
        It is also called when a string menu item is selected.

        Args:
            messageId (int): The ID of the gadget that triggered the event.
            bc (c4d.BaseContainer): The original message container.

        Returns:
            bool: False if there was an error, otherwise True.
        """
        # User click on Creat
        if messageId == self.ID_BUTTON_CREATE:
            #check path and creat rdkit mol
            print(self.GetFilename(self.ID_FILEPATH), "Create MOL")
            mol = self.ReadMol(self.GetFilename(self.ID_FILEPATH))
            if mol == False:
                return True
            # print(self.GetFilename(self.ID_FILEPATH), "Create MOL")
            geo, geo_bond = self.BuildPDTGeo(mol)
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
            
            pos = geo.GetP()
            pmesh = C4DFunction.postomesh(pos)
            pmesh[c4d.ID_BASELIST_NAME] = "PositionMesh"
            pmesh.InsertUnder(rootobj)
            t = c4d.BaseTag(c4d.Tpython)

            # bc_buildtype = c4d.GetCustomDataTypeDefault(c4d.DTYPE_LONG)
            bc_buildtype = c4d.gui.QuickTabCustomGui()

            # bc_buildtype[c4d.DESC_NAME] = "Build type"
            # bc_buildtype[c4d.DESC_CUSTOMGUI] = c4d.CUSTOMGUI_QUICKTAB
            bc_buildtype.AppendString(1, "a")
            bc_buildtype.AppendString(2, "b")
            bc_buildtype[c4d.QUICKTAB_BAR] = False
            descId_buildtype = t.AddUserData(bc_buildtype)
            t[descId_buildtype] = 1
            #iterate atom symbols and build hierarchy
            indnum = 0
            rad = 1

            for i in geo.pointgroups:
                # create userdata on python tag-----------------------------------------------------------
                # use template to set type and name
                bc_radius = c4d.GetCustomDataTypeDefault(c4d.DTYPE_REAL)
                bc_radius[c4d.DESC_NAME] = "Radius " + i
                bc_radius[c4d.DESC_MIN] = 0
                bc_radius[c4d.DESC_MAX] = 200
                # get id
                descId_radius = t.AddUserData(bc_radius)
                # set value
                t[descId_radius] = parameters[2*indnum]*rad

                bc_color = c4d.GetCustomDataTypeDefault(c4d.DTYPE_COLOR)
                bc_color[c4d.DESC_NAME] = "Color " + i
                descId_color = t.AddUserData(bc_color)
                t[descId_color] = parameters[2*indnum+1]

                #creat atoms---------------------------------------------------------------------------------
                #creat point selection tags
                seltag = c4d.BaseTag(c4d.Tpointselection)
                seltag[c4d.ID_BASELIST_NAME] = i

                #the returned BaseSelect obj can be modified and automatically updated back to selection tag
                sel = seltag.GetBaseSelect()
                #select points from mark list, 1 is selected, 0 is unselected
                marklist = geo.GetPointGroup(i)
                sel.SetAll(marklist)
                
                #attach tag to mesh obj
                pmesh.InsertTag(seltag)

                #setup cloner
                sph = c4d.BaseObject(5160)
                sph[c4d.ID_BASELIST_NAME] = i
                sph[c4d.PRIM_SPHERE_RAD] = 1
                sph[c4d.PRIM_SPHERE_TYPE] = 4
                sph[c4d.PRIM_SPHERE_SUB] = 10
                sphphone = c4d.BaseTag(5612)
                sphphone[c4d.PHONGTAG_PHONG_ANGLELIMIT] = 1
                sph.InsertTag(sphphone)
                cln = c4d.BaseObject(1018544)
                cln[c4d.ID_BASELIST_NAME] = i
                cln[c4d.MGCLONER_VOLUMEINSTANCES_MODE] = 1
                cln[c4d.ID_MG_TRANSFORM_SCALE] = c4d.Vector(t[descId_radius])
                cln[c4d.ID_MG_TRANSFORM_COLOR] = t[descId_color]
                cln[c4d.ID_MG_MOTIONGENERATOR_MODE] = 0
                cln[c4d.MG_OBJECT_LINK] = pmesh
                cln[c4d.MG_POLY_MODE_] = 0
                cln[c4d.MG_POLY_SELECTION] = i

                sph.InsertUnder(cln)
                cln.InsertUnderLast(rootatom)

                #creat bond splines------------------------------------------------------------------------
                #iterate bonds and get position list for atom type i
                #creat empty spline
                spl = c4d.BaseObject(5101)
                linepos = []
                for prim in geo_bond.primitives:
                    if prim.symbol == i:
                        for vertex in prim.vertices:
                            linepos.append(geo_bond.GetP()[vertex.pointnumber])
                #set point count and segment count
                spl.ResizeObject(len(linepos), int(len(linepos)/2))

                #set point positions and assign point to segment
                spl.SetAllPoints(linepos)
                for j in range(int(len(linepos)/2)):
                    spl.SetSegment(j, 2, False)
                
                #setup sweep
                spl[c4d.ID_BASELIST_NAME] = i
                spl.Message(c4d.MSG_UPDATE)
                circle = c4d.BaseObject(5181)
                circle[c4d.PRIM_CIRCLE_RADIUS] = 10
                sweep = c4d.BaseObject(5118)
                sweepphone = c4d.BaseTag(5612)
                sweepphone[c4d.PHONGTAG_PHONG_ANGLELIMIT] = 1
                sweep.InsertTag(sweepphone)
                sweep[c4d.ID_BASEOBJECT_USECOLOR] = 2
                sweep[c4d.ID_BASEOBJECT_COLOR] = t[descId_color]
                sweep[c4d.ID_BASELIST_NAME] = i
                spl.InsertUnder(sweep)
                circle.InsertUnder(sweep)
                sweep.InsertUnderLast(rootbond)

                indnum = indnum +1
            rootobj.InsertTag(t)
            doc.InsertObject(rootobj)

            t[c4d.TPYTHON_CODE] = '''
from typing import Optional
import c4d

doc: c4d.documents.BaseDocument # The document evaluating this tag
op: c4d.BaseTag # The Python scripting tag
flags: int # c4d.EXECUTIONFLAGS
priority: int # c4d.EXECUTIONPRIORITY
tp: Optional[c4d.modules.thinkingparticles.TP_MasterSystem] # Particle system
bt: Optional[c4d.threading.BaseThread] # The thread executing this tag

def main() -> None:
    # Called when the tag is executed. It can be called multiple time per frame. Similar to TagData.Execute.
    # Write your code here
    root = op.GetObject()
    root_bond = root.GetDownLast().GetPred().GetChildren()
    root_atom = root.GetDownLast().GetChildren()
    # print(root_bond,root_atom)
    t = root.GetTag(c4d.Tpython)
    for i in range(len(root_atom)):
        root_bond[i][c4d.ID_BASEOBJECT_COLOR] = t[c4d.ID_USERDATA, 2*i+2]
        root_atom[i][c4d.ID_MG_TRANSFORM_SCALE] = c4d.Vector(t[c4d.ID_USERDATA, 2*i+1])
        root_atom[i][c4d.ID_MG_TRANSFORM_COLOR] = t[c4d.ID_USERDATA, 2*i+2]

'''
            c4d.EventAdd()
            return True
        return True
    
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

    

