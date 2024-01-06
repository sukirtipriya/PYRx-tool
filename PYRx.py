
import os, sys, subprocess, pickle
from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4LigandPreparation, AD4ReceptorPreparation, AD4FlexibleReceptorPreparation
from AutoDockTools.GridParameters import GridParameter4FileMaker
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
from MolKit.pdbWriter import PdbWriter
from MolKit.pdbParser import PdbqtParser
from MolKit.molecule import BondSet
try:
    from enthought.preferences.ui.api import PreferencesPage
    from enthought.traits.api import Range, Int, Str, File, Directory, Trait, List, Bool, HTML
    from enthought.preferences.api import get_default_preferences
    from enthought.traits.ui.api import View, Item, RangeEditor, Group, CheckListEditor, HGroup, HTMLEditor
except:
    from apptools.preferences.ui.api import PreferencesPage
    from traits.api import Range, Int, Str, File, Directory, Trait, List, Bool, HTML
    from apptools.preferences.api import get_default_preferences
    from traitsui.api import View, Item, RangeEditor, Group, CheckListEditor, HGroup, HTMLEditor    
import wx
from utils import rcFolder, which
from miscTraits import PositiveInt
from MolKit.protein import Protein, Chain, Residue
from MolKit.molecule import Atom
import shutil, copy

class AutoDockPreferencesPage(PreferencesPage):
    workspace = Directory()

    autodock = File('autodock4', filter=["*"])
    autogrid = File('autogrid4', filter=["*"])
    vina = File('vina', filter=["*"])
    cpu_num = PositiveInt()
    max_time = PositiveInt(1000000)
    URI = Str()
    executionMode = Range( 0, 3 ) 

    #Used to limit the number of docking to add to in AnalyzeVinaPage 
    limit_docking_for_gui_analysis = Bool(False)    
    number_of_top_results_to_display = Range( 0, sys.maxint, 1000 ) #depens of available memory
        
    
    #### Traits UI views ######################################################

    view = View(Group(
                        Group('autodock', 'autogrid', 'vina',
                        Item(name='workspace'),
                        Item(name='cpu_num',  label="Available CPUs:"),
                        Item(name='max_time',  label="Max Runtime (m):"),),
                        Item(name='limit_docking_for_gui_analysis', label="Limit Docking for GUI Analysis"),
                        Item(name='number_of_top_results_to_display', visible_when ='object.limit_docking_for_gui_analysis==True'), 
                                   label="AutoDock Preferences"
                    ),
                )
    
    def _autodock_changed(self, new):
        startPage = wx.GetApp().GetTopWindow().autodockWiz.book.GetPage(0)
        if new and which(new) and which(self.autodock):
            startPage.rb.EnableItem(0, True)
        else:
            startPage.rb.EnableItem(0, False)            
        startPage.Update()
        startPage.rb.SetString(0,which(self.autodock))
        
    def _autogrid_changed(self, new):
        startPage = wx.GetApp().GetTopWindow().autodockWiz.book.GetPage(0)
        if new and which(new) and which(self.autodock):
            startPage.rb.EnableItem(0, True)
        else:
            startPage.rb.EnableItem(0, False)
        startPage.Update()
        
    def _vina_changed(self, new):
        startPage = wx.GetApp().GetTopWindow().vinaWiz.book.GetPage(0)
        if new and which(new) and which(self.vina):
            startPage.rb.EnableItem(0, True)
        else:
            startPage.rb.EnableItem(0, False)
        startPage.rb.SetString(0,which(self.vina))
        startPage.Update()     
        
    def _workspace_changed(self, new):
        if new == self.workspace: return
        frame = wx.GetApp().GetTopWindow()  
        frame.SetAllCursors(wx.StockCursor(wx.CURSOR_WAIT))
        frame.Refresh()  
        frame.navigator.SetSelection(frame.navigator.GetPageIndex(frame.autodockNav.autodockTree))
        frame.vsModel.__init__()
        frame.autodockNav.autodockTree.ligandTree.tree.DeleteAllItems()
        frame.autodockNav.autodockTree.ligandTree.BuildTree(frame.vsModel.ligandsFolder)
        frame.autodockNav.autodockTree.macromoleculeTree.tree.DeleteAllItems()
        frame.autodockNav.autodockTree.macromoleculeTree.BuildTree(frame.vsModel.macromoleculesFolder)
        frame.SetAllCursors(wx.NullCursor)

autodockPreferencesPage = AutoDockPreferencesPage()

class LigandPreparationPage(PreferencesPage):
    category = 'AutoDock'
    repairs = List( editor = CheckListEditor( 
                           values = [ 'bonds', 'hydrogens',], 
                           cols   = 2 ) ) 
    inactivate_all_torsions = Bool()
    limit_torsions = Bool(False)    
    limit_torsion_number = Range( 0, 120, 1 )
    #### Traits UI views ######################################################

    view = View(Group(  Item(name='repairs', label = "Type(s) of repairs to make", style='custom'),
                        Item(name='inactivate_all_torsions',),             
                        Item(name='limit_torsions', ),
                        Item(name='limit_torsion_number', label="Number of torsions", visible_when ='object.limit_torsions ==True'), 
                        label="Ligand Preparation Preferences"
                    ),

                )           
ligandPreparationPage = LigandPreparationPage()

class ReceptorPreparationPage(PreferencesPage):
    category = 'AutoDock'

    cleanup = List( editor = CheckListEditor( 
                                             values = [ 'nphs', 'lps', 'waters', 'nonstdres'], cols = 2,
                                             ), value = ['nphs', 'lps', 'waters'], 
                   ) 
    txt = Str("""
 Nphs - merge charges and remove non-polar hydrogens
 Lps - merge charges and remove lone pairs
 Waters - remove water molecules 
 Nonstdres - remove chains composed entirely of residues of types 
                  other than the standard 20 amino acids""")
    #### Traits UI views ######################################################
    view = View(Group(  Item(name='cleanup', label = "Type(s) of changes to make", style='custom'),
                        label="Receptor Preparation Preferences"
                     ),
                Item(name="txt", style='readonly', show_label=False)
                )           
receptorPreparationPage = ReceptorPreparationPage()

class AutoDockHPCPreferencesPage(PreferencesPage):
    category = 'AutoDock'

    key_pair = File('key.pem', filter=["*"])
    host_name = Str()
    help = HTML("Click on the link below for help on how to set this up. <br><br><a href='http://goo.gl/oMSzWu'>http://goo.gl/oMSzWu</a>")

    #### Traits UI views ######################################################
    
    view = View(Group(
                            Item(name='key_pair'),
                            Item(name='host_name', label = "Host Name or IP Address"),
                            Item(name='help', editor=HTMLEditor(open_externally=True), style='readonly' , height=100),
                            label="HPC Cluster Preferences"
                            ),
                        )
    
hpc_modal_view = View(Group(
                        Item(name='key_pair'),
                        Item(name='host_name', label = "Host Name or IP Address"),
                        Item(name='help', editor=HTMLEditor(open_externally=True), style='readonly' , height=80),
                        label="HPC Cluster Preferences"
                        ),
                    buttons = ['OK', 'Cancel' ], kind = 'livemodal', width=500                   
                    )    
autodockHPCPreferencesPage = AutoDockHPCPreferencesPage()

class AutoDockInverseVSPreferencesPage(PreferencesPage):
    category = 'AutoDock'
    max_macroMols = PositiveInt(100) #arbitrary number
    text = Str("""For Inverse Virtual Screening with Vina, load 
maximum this number of Macromolecules. Grid dimensions 
for remaining Macromolecules are based on the last 
loaded Macromolecule. This number can be changed by  
Users depending on available memory (RAM) and the size 
of Macromolecules. On a workstation with 8GB of RAM, 
for instance, you can load around 100 Macromolecules 
each with 2500 atoms before it runs out of memory.""")
    #### Traits UI views ######################################################
    view = View(Group(
                            Item(name='max_macroMols', label = "Max Macromolecules"),
                            label="Inverse Virtual Screening Preferences"
                            ),
                Item(name='text', label="Help",  style = 'readonly'),
                        )           
autodockInverseVSPreferencesPage = AutoDockInverseVSPreferencesPage()

class VSModel:
    """
    Used for setting up and running Virtual Screening.
Example:
>>> vs = VSModel(ligands=[ind.pdb],macroMolecule='hsg1.pdb')
#creates user.home/AutoDock4Data/
#                               Macromolecules (stores Macromolecules)
#                               Ligands        (stores Ligand)
>>> vs.run()
#runs AutoGrid followed by AutoDock and stores results in
user.home/AutoDock4Data/
                        hsg1/
                            hsg1.A.map
                            ...
                            hsg1.gpf    (stores AutoGrid runs)
                            hsg1.pdbqt
                            ind/
                                ind.dlf  (stores AutoDock Runs)
                                ind.dpf
user.home/AutoDock4Data/                                        
                       ind.pdbqt
    """

    def __init__(self, ligandPaths=[], macromoleculePath=None):
        """
        basePath - path for storing the data, if None user.home/AutoDock4Data is used.
self.ligandPaths - lists of ligand paths.
self.macromoleculePath  - path to Macromolecule.
        """
        pref = get_default_preferences()
        basePath = pref.get('AutoDock.workspace')
            
                
        #Note Folder creation can be done elsewhere
        if not os.path.isdir(basePath):
            os.mkdir(basePath)
        self.macromoleculesFolder = os.path.join(basePath,'Macromolecules')
        if not os.path.isdir(self.macromoleculesFolder):
            os.mkdir(self.macromoleculesFolder)            
        self.ligandsFolder = os.path.join(basePath,'Ligands')
        if not os.path.isdir(self.ligandsFolder):
            os.mkdir(self.ligandsFolder) 
        self.etcFolder = os.path.join(basePath,'etc')
        if not os.path.isdir(self.etcFolder):
            os.mkdir(self.etcFolder) 
        
        self.basePath = basePath        
        self.ligandPaths = ligandPaths
        self.macromoleculePath = macromoleculePath
        self.LPO_list = []        
        self.ligand_types = set()
        
    def CheckMaps(self):
        "Check to see if all maps dimensions are the same"  
        ligandTypes = list(self.ligand_types)
        firstMap = os.path.join(self.macromolecule.receptorFolder, self.macromolecule.receptor_stem+'.'+ligandTypes[0]+'.map')
        firstMapFile = open(firstMap)
        lineList = [] 
        numberCheck = 6
        for i in range(numberCheck):
            lineList.append(firstMapFile.readline())
        for ligandType in ligandTypes[1:]:
            filePath = os.path.join(self.macromolecule.receptorFolder, self.macromolecule.receptor_stem+'.'+ligandType+'.map')
            mapFile = open(filePath)
            for i in range(numberCheck):
                line = mapFile.readline()
                if line != lineList[i]:
                    mapFile.close()
                    firstMapFile.close()
                    return False      
            mapFile.close()
        firstMapFile.close()
        return True
     
    def Run(self):
        self.PrepareAllLigands()
        self.PrepareReceptor()
        self.PrepareGPF()
        #self.runAutoGrid()
        #os.waitpid(self.AutoGridProcess.pid,0)
        self.RunAllDocking()
        
    def RrepareAllLigands(self):
        for ligandPath in self.ligandsPaths:
            if not os.path.exists(ligandPath):
                print "ligand -"  + ligandPath + " - does not exists."
            else:
                self.PrepareLigandFile(ligandPath)
        self.GetMolDict()
        
    def PrepareLigandFile(self, ligandFile):
        mols = Read(ligandFile)
        if len(mols)>1:
            print "%d molecules in %s"%(len(mols), ligandFile)
            print "%s will use the first molecule as ligand"%(self.__module__)
        mol = mols[0]
        mol.buildBondsByDistance()
        outputfilename = os.path.join(self.ligandsFolder, mol.name + ".pdbqt")
        LPO = AD4LigandPreparation(mol, outputfilename=outputfilename)
        self.LPO_list.append(LPO)
        self.PrepareLigandMol(mol)

    def PrepareLigandMol(self, mol, charges_to_add='gasteiger', use_=False):
        if ligandPreparationPage.limit_torsions:
            limit_torsions =  ligandPreparationPage.limit_torsion_number
        else:
            limit_torsions = False
        mol.name = mol.name.replace('*','_')
        outputfilename = os.path.join(self.ligandsFolder, mol.name + ".pdbqt")
        if use_:
            same_outCounter = 1 #this is in case outputfilename exists
            while os.path.exists(outputfilename):
                outputfilename = os.path.join(self.ligandsFolder, mol.name + "_"+ str(same_outCounter)+ ".pdbqt")
                same_outCounter += 1
                
        LPO = AD4LigandPreparation(mol, outputfilename=outputfilename, charges_to_add=charges_to_add,
                                   repairs = '_'.join(ligandPreparationPage.repairs),
                                   inactivate_all_torsions = ligandPreparationPage.inactivate_all_torsions,
                                   limit_torsions = limit_torsions)
#        conectRecords = None
#        if hasattr(mol, '_openBebel'):
#            writer = PdbWriter()
#            writer.defineConnectSection(mol.allAtoms, bondOrigin='File')
#            conectRecords = writer.recordsToWrite['CONECT']
#            
        pickleFileName = os.path.join(self.etcFolder, mol.name + ".pkl")
        #This is needed in order to store autodock_element or TORSDOF in a file that we can access without reading the hole file
        molDict = {'autodock_element':set(LPO.molecule.allAtoms.autodock_element)}
        molDict['TORSDOF']  = LPO.molecule.TORSDOF
        LPO.molecule.getCenter()
        molDict['center'] = LPO.molecule.center
#        if conectRecords:
#            molDict['CONECT'] = conectRecords
        pickle.dump(molDict, open(pickleFileName, 'w'))
        if len(LPO.molecule.allAtoms) > len(LPO.writer.writtenAtoms):
            txt = "AD4LigandPreparation wrote less atoms that present in the molecule: "+mol.name
            txt += "\nAD4LigandPreparation wrote less atoms that present in the molecule: "+mol.name
            txt += "\nLocation:"+__file__ + ":PrepareLigandMol"
            try:
                 wx.GetTopLevelWindows()[0].log.warning(txt)
            except:
                print txt
            
            mol = Read(outputfilename)[0]
            return mol, molDict
        return LPO.molecule,  molDict

    def CreateMolDict(self):
        """Prepares molecule dictionary with 
        key   - molecules name
        value - dictionary ({'autodock_element':...) taken from a pickle file generated with PrepareLigandMol
        """
        molDict = {}
        ligand_types = set()

        for ligand in self.ligands:
            name = os.path.splitext(os.path.split(ligand)[1])[0]
            pickleFileName = os.path.join(self.etcFolder, name + ".pkl")    
            if os.path.exists(pickleFileName):
                dict = pickle.load(open(pickleFileName))
            else:
                try:
                    parser = PdbqtParser(ligand)
                    mol = parser.parse()[0]
                    pickleFileName = os.path.join(self.etcFolder, mol.name + ".pkl")
                    #This is needed in order to store autodock_element or TORSDOF in a file that we can access without reading the hole file
                    dict = {'autodock_element':set(mol.allAtoms.autodock_element)}
                    dict['TORSDOF']  = mol.TORSDOF
                    mol.getCenter()
                    dict['center'] = mol.center
                    pickle.dump(dict, open(pickleFileName, 'w'))
                except Exception,e:
                    wx.GetApp().GetTopWindow().log.warn("Error parsing "+ligand+"\n"+str(e))
                    continue
            dict['name'] = name.encode()
            molDict[name] = dict
            ligand_types = ligand_types.union(dict['autodock_element'])
        self.molDict = molDict
        self.ligand_types = ligand_types


    def PrepareReceptorMol(self, mol):        
        receptorFolder = os.path.join(self.macromoleculesFolder, mol.name)
        if not os.path.isdir(receptorFolder):
            os.mkdir(receptorFolder)
        self.receptorFolder = receptorFolder
        outputfilename = os.path.join(receptorFolder, mol.name + ".pdbqt")
        RPO = AD4ReceptorPreparation(mol, outputfilename=outputfilename, 
                                     cleanup='_'.join(receptorPreparationPage.cleanup))
        self.macromolecule = RPO.molecule
        self.macromolecule.allAtoms = RPO.molecule.chains.residues.atoms #this is needed to update allAtoms
        self.macromoleculePath = outputfilename
        return RPO
    
    def PrepareFlexReceptor(self, flexRes):
        "Prepares _rigid.pdbqt and _flex.pdbqt files"        
        name = flexRes[0].top.name
        receptorFolder = os.path.join(self.macromoleculesFolder, name+"_flex")
        if not os.path.isdir(receptorFolder):
            os.mkdir(receptorFolder)
        self.receptorFolder = receptorFolder               
        rigid_filename = os.path.join(receptorFolder, name + "_rigid.pdbqt")
        AD4ReceptorPreparation(flexRes[0].top, outputfilename=rigid_filename, cleanup='_'.join(receptorPreparationPage.cleanup))
        flexres_filename = os.path.join(receptorFolder, name+"_flex.pdbqt")
        flex_res_lines = []
        for residue in flexRes: #handle one residue at a time. 
            prot = Protein()
            allAtoms = copy.copy(residue.atoms.copy())
            chain = Chain(residue.parent.name,  prot, top=prot)
            res = Residue(residue.type, number=residue.number, parent=chain, top=prot)
#            res.atoms = residue.atoms
            for atom in allAtoms:
                res.adopt(atom)
                for b in atom.bonds:
                    if not b.atom1 in allAtoms or not b.atom2 in allAtoms:
                        atom.bonds.remove(b)
            prot.allAtoms = allAtoms
            prot.allAtoms.top = prot
            prot.parser = flexRes[0].top.parser
            prot.levels = [Protein, Chain, Residue, Atom]
            AD4FlexibleReceptorPreparation(prot, residues=[res], rigid_filename=rigid_filename+"_", flexres_filename=flexres_filename)
            flex_res_lines.extend(open(flexres_filename).readlines())
        os.remove(rigid_filename+"_")
        flexres_file = open(flexres_filename, 'w')
        for line in flex_res_lines:
            flexres_file.write(line)
        flexres_file.close()
        cmp_flex_res_lines = []
        for line in flex_res_lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                cmp_flex_res_lines.append(line[15:])
        rigid_file_lines = open(rigid_filename).readlines()
        rigid_file = open(rigid_filename, 'w')
        for line in rigid_file_lines:
            if not line[15:] in cmp_flex_res_lines:
                rigid_file.write(line)
        rigid_file.close()
        self.macromoleculePath = rigid_filename
        self.flexres_filename = flexres_filename
    
    def PrepareGPF(self, gridParameters=None):
        gpfm = GridParameter4FileMaker()
        gpfm.receptor = self.macromolecule
        #set_receptor part
        receptor_types = gpfm.getTypes(gpfm.receptor)
        ad4_typer = AutoDock4_AtomTyper()
        ad4_typer.setAutoDockElements(gpfm.receptor)
        receptor_types = set(gpfm.receptor.allAtoms.autodock_element)
        #gpo.set_receptor4 part

        gpfm.gpo['receptor_types']['value'] = ' '.join(list(receptor_types))
        gpfm.gpo['ligand_types']['value'] = ' '.join(list(self.ligand_types))
        gpfm.gpo.receptor_filename = self.macromolecule.receptor_filename
        receptor_stem = self.macromolecule.receptor_stem
        gpfm.gpo.receptor_stem = receptor_stem
        gpfm.gpo['receptor']['value'] = self.macromolecule.receptor_filename
        gpfm.gpo['gridfld']['value'] = receptor_stem + '.maps.fld'
        gpfm.gpo['elecmap']['value'] = receptor_stem + '.e.map'
        gpfm.gpo['dsolvmap']['value'] = receptor_stem + '.d.map'
        if gridParameters:
            for key in gridParameters:
                kw = {key:gridParameters[key]}
                gpfm.set_grid_parameters(**kw)
        inputFile = receptor_stem + ".gpf"
        output_gpf_filename = os.path.join(self.macromolecule.receptorFolder, inputFile)
        self.gpf_filename = output_gpf_filename
        self.gpfm = gpfm
        gpfm.write_gpf(output_gpf_filename)
        outputFile = self.macromolecule.name + ".glg"
        self.gridCommand = [autodockPreferencesPage.autogrid, "-p", inputFile, "-l", outputFile]
        self.glgOutput = os.path.join(self.macromolecule.receptorFolder, outputFile)
        
    def RunAutoGrid(self):
        cwd = os.path.split(self.macromoleculePath)[0]
        outputFile = self.macromolecule.name + ".glg"
        cmd = autodockPreferencesPage.autogrid + " -p " + self.gpf_filename + " -l " + outputFile
        self.AutoGridProcess = subprocess.Popen(cmd, stdin=subprocess.PIPE, 
                                                stdout=subprocess.PIPE, 
                                                stderr=subprocess.PIPE,
                                                cwd=cwd, shell=True)
        
    def RunAllDocking(self):
        for LPO in self.LPO_list:
            self.prepareDPF(LPO)
            self.runAutoDock()
            
    def QsubAllDocking(self):
        for LPO in self.LPO_list:
            self.prepareDPF(LPO)
            self.qsubAutoDock()
        
    def PrepareDPF(self, ligandParameters, docking_algorithm_parameter_list):
        basename = ligandParameters['name']+'.pdbqt'
        self.ligand = basename
        basename = os.pardir + os.sep + os.pardir + os.sep + os.pardir + os.sep + "Ligands" + os.sep + basename
        self.dpo['rmsref']['value'] = basename
        self.dpo['move']['value'] = basename        
        self.dpo['torsdof4']['value'][0] = ligandParameters['TORSDOF']
        self.dpo['ndihe']['value'] = ligandParameters['TORSDOF']
        ADelement = list(ligandParameters['autodock_element'])
        self.dpo['ligand_types']['value'] = ' '.join(list(self.ligand_types))
        dockingFolder = os.path.join(self.macromolecule.receptorFolder, ligandParameters['name'])
        if not os.path.isdir(dockingFolder):
            os.mkdir(dockingFolder)        
        dockingFolder = os.path.join(self.macromolecule.receptorFolder, ligandParameters['name'])
        if self.macromolecule.receptor_stem.endswith('_rigid'):
            self.dpo['flexres']['value'] = os.pardir+os.path.sep+self.macromolecule.receptor_stem.replace('_rigid','_flex')+".pdbqt"
            self.dpo['flexres_flag']['value'] = True            
        cen = ligandParameters['center']
        self.dpo['about']['value'] =  [round(cen[0],4), round(cen[1],4),\
                                        round(cen[2],4)]
        dpf_filename = self.macromolecule.receptor_stem + '_' + ligandParameters['name'] + '.dpf'
        dpf_filename = os.path.join(dockingFolder, dpf_filename)
        self.dpo.set_receptor(self.macromolecule.receptor_filename)
        self.dpo.write42(dpf_filename, docking_algorithm_parameter_list)
        txt = open(dpf_filename).read()
        txt = txt.replace(self.macromolecule.receptor_stem+".", os.pardir+os.sep+self.macromolecule.receptor_stem+".")
        txt = txt.replace(basename, basename+" ")  #TODO: remove this after PyRx 1.0 release
        
        open(dpf_filename,'w').write(txt)
        self.dockingFolder = dockingFolder
        self.dpf_filename = os.path.basename(dpf_filename)
        self.full_dpf_filename = dpf_filename
        outputFile = self.dpf_filename.replace('.dpf','.dlg')
        self.dockCommand = [autodockPreferencesPage.autodock, "-p", self.dpf_filename, "-l", outputFile]
        self.dlgOutput = os.path.join(dockingFolder, outputFile)
        
    def RunAutoDock(self):
        cmd = autodockPreferencesPage.autodock + " -p " + self.dpf_filename + " -l " + self.dpf_filename.replace('.dpf','.dlg')
        if sys.platform != 'win32':
            cmd = 'ulimit -s unlimited;' + cmd
        AutoDockProcess = subprocess.Popen(cmd, stdin=subprocess.PIPE, 
                                                stdout=subprocess.PIPE, 
                                                stderr=subprocess.PIPE,
                                                cwd=self.dockingFolder, shell=True)
        return AutoDockProcess
    
    def QsubAutoDock(self):
        txt = "ulimit -s unlimited\n"
        txt += "cd " + self.dockingFolder + "\n"
        txt += autodockPreferencesPage.autodock + " -p " + self.dpf_filename + " -l " + self.dpf_filename.replace('.dpf','.dlg')
        jobFile = self.dpf_filename.replace('.dpf','_AD')
        open(jobFile,'w').write(txt)
        cmd = "chmod +x " + jobFile + "\n"
        cmd += "qsub -l cput=23:00:00 -l nodes=1:ppn=1 -l walltime=23:30:00 -l mem=512mb " + jobFile
        #jobIDsName = jobFile + ".jobIDs"
        #cmd += " >>" + jobIDsName
        os.system(cmd)
