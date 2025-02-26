global cd, Lx, Ly, T, nOrder_contact, nOrder_potts, Time, MCS0
global LamV, tV
global pP, pD, pC
global E_pd, E_dd, E_pp, E_ll, E_lp, E_ld, E_ml, E_mp, E_md, E_mm, E_mc, E_mt, E_tt, E_cc, E_cp, E_cd, E_cl
global J_pd, J_other

#PARAMETERS:
#volume of each cell=144 for square lattice and volume=169 for hexagonal lattice.
cd=12           #Typical cell diameter (cd=13 for hexagonal lattice)
nx = 4          #Number of cells along x-direction
ny = 4          #Number of cells along y-direction
Ly = nx*cd      #for periodic boundary 
Lx = nx*cd      #for periodic boundary 
#Ly = nx*cd           #for left boundary signal configuration 
#Lx = nx*cd + 20      #for left boundary signal configuration 
#
#POTTS PARAMETERS:
T=15            #Temperature (T=20 for hexagonal lattice)
nOrder_contact=4 #Distance of interaction (nOrder=3 for hexagonal lattice)
nOrder_potts=3  #Distance of interaction (nOrder=2 for hexagonal lattice)
Time=100000     #Total number of MCS
#
#VOLUME/SURFACE PARAMETERS:
LamV=12         #Lambda Volume
tV=144          #Target Volume (tV=169 for hexagonal lattice)
#
#SUBCELLULAR PARAMETERS:
#volume fractions 
pP=0.15        #Proximal domain (for square lattice)
pD=0.15        #Distal domain (for square lattice)
pC=0.5         #Cytoplamsic domain (for square lattice)
# pP=0.142     # Proximal domain (for hexagonal lattice)
# pD=0.142     #Distal domain (for hexagonal lattice)
# pC=0.574     # Cytoplamsic domain (for hexagonal lattice)
#
#EXTERNAL CONTACT ENERGIES:
E_pd = 7
E_pp = 12
E_dd = E_pp      
E_lp = 12
E_ld = E_lp
E_ll = 14
E_ml = E_ll
E_mp = 10
E_md = E_mp
E_mm = 0
E_mc = 100
E_mt = 5
E_tt = E_mt
E_cc = 50
E_cp = 30
E_cd = E_cp
E_cl = E_cp
#
#INTERNAL CONTACT ENERGIES         #Internal contact energy between proximal and distal
J_pd = 30
J_other = 1
#
t_relax_growth = 10000             #time taken by the cells to grow and relax their boundaries after placing the initial seeds
#
#CELL PROLIFERATION PARAMETERS:
growth_rate= 0.001                 #rate of growth of cells during proliferation
t_relax = 10000                    #relaxation time for cell proliferation (set to 100,00,000 for uniform cell proliferation)                 
#
MCS0 = 0                           #set to zero when simulation should restart from a .piff file (used for resuimg simulations from .piff files)
#
def configure_simulation():

    from cc3d.core.XMLUtils import ElementCC3D
    
    CompuCell3DElmnt=ElementCC3D("CompuCell3D",{"Revision":"0","Version":"4.3.1"})
    
    MetadataElmnt=CompuCell3DElmnt.ElementCC3D("Metadata")
    
    MetadataElmnt.ElementCC3D("NumberOfProcessors",{},1)
    MetadataElmnt.ElementCC3D("DebugOutputFrequency",{},"0")
    
    PottsElmnt=CompuCell3DElmnt.ElementCC3D("Potts")
    
    PottsElmnt.ElementCC3D("Dimensions",{"x":Lx,"y":Ly,"z":"1"})
    PottsElmnt.ElementCC3D("Steps",{},Time)
    PottsElmnt.ElementCC3D("Temperature",{},int(T))
    PottsElmnt.ElementCC3D("NeighborOrder",{},nOrder_potts)
    PottsElmnt.ElementCC3D("Boundary_x",{},"Periodic")          #comment this line when using left boundary orienting signal
    PottsElmnt.ElementCC3D("Boundary_y",{},"Periodic")
    PottsElmnt.ElementCC3D("LatticeType",{},"Square")
    
    PluginElmnt=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CellType"})
    
    # Listing all cell types in the simulation
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"0","TypeName":"Medium"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"1","TypeName":"cyto"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"2","TypeName":"prox2"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"3","TypeName":"dist2"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"4","TypeName":"lat"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"5","TypeName":"test"})        #the seeds of initial cells are defined as type "test"
    
    CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"BoundaryPixelTracker"})
    
    CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"VolumeLocalFlex"})
    
    PluginElmnt_1=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CenterOfMass"})
    
    PluginElmnt_2=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"NeighborTracker"})
    
    PluginElmnt_4=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Contact"})
    # Specification of external contact energies
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"Medium"},E_mm)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"cyto"},E_mc)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"prox2"},E_mp)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"dist2"},E_md)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"lat"},E_ml)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"test"},E_mt)
    
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"test","Type2":"test"},E_tt)
    
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"cyto","Type2":"cyto"},E_cc)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"cyto","Type2":"prox2"},E_cp)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"cyto","Type2":"dist2"},E_cd)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"cyto","Type2":"lat"},E_cl)
    
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"prox2","Type2":"prox2"},E_pp)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"prox2","Type2":"dist2"},E_pd)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"prox2","Type2":"lat"},E_lp)
    
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"dist2","Type2":"dist2"},E_dd)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"dist2","Type2":"lat"},E_ld)
    
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"lat","Type2":"lat"},E_ll)
    PluginElmnt_4.ElementCC3D("NeighborOrder",{},nOrder_contact)

    
    PluginElmnt_5=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"ContactInternal"})
    # Specification of internal contact energies
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"cyto","Type2":"prox2"},J_other)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"cyto","Type2":"test"},0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"cyto","Type2":"dist2"},J_other)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"cyto","Type2":"lat"},J_other)
    
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"prox2","Type2":"dist2"},J_pd)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"prox2","Type2":"lat"},J_other)

    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"dist2","Type2":"lat"},J_other)
    PluginElmnt_5.ElementCC3D("NeighborOrder",{},nOrder_contact)
    
    # SteppableElmnt=CompuCell3DElmnt.ElementCC3D("Steppable",{"Type":"PIFInitializer"})        #used to resume simulations from .piff files
    # SteppableElmnt.ElementCC3D("PIFName",{},"PiffFile_100000.piff")
    
    CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)
    
from cc3d import CompuCellSetup
        

configure_simulation()            

            

from PCPSteppables import *

CompuCellSetup.register_steppable(steppable=InitialConditionsSteppable(10,cd,tV,LamV,pP,pD,pC))                #Initial conditions used for cell proliferation
# CompuCellSetup.register_steppable(steppable=periodic_boundary(10,cd,tV,LamV,pP,pD,pC,nx,ny,t_relax_growth))     Initial conditions used for periodic boundary conditions
# CompuCellSetup.register_steppable(steppable=left_boundary(10,cd,tV,LamV,pP,pD,pC,nx,ny,t_relax_growth))        Initial conditions used for left boundary orienting signal configuration
CompuCellSetup.register_steppable(steppable=Volume_steppable(10,tV,LamV))                                       #Defines the volume constraints of the cells
CompuCellSetup.register_steppable(steppable=CellPol(10,cd))                                                     #Computes and assigns polarity vectors to CYTO cells based on the relative positions of proximal and distal cells within a cluster.
CompuCellSetup.register_steppable(steppable=GraphGlobalPol(10,cd,MCS0))                                         #Computes the global polarization (Phi) by averaging individual cell polarity vectors and saves the data into a .txt file
CompuCellSetup.register_steppable(steppable=OrientationField(10))                                              #Visualizes cell polarity by creating vector fields
#CompuCellSetup.register_steppable(steppable=VectorField(1000))                                                 #Saves vector field data into a .txt file
#CompuCellSetup.register_steppable(steppable=GraphPol_x_Dist(1000,cd))                                          #Calculates distance-dependent polarization for periodic boundary configurations
CompuCellSetup.register_steppable(steppable=GraphPol_x_Dist_left_boundary(1000,cd,t_relax_growth))               #Calculates distance-dependent (with respect to the left boundary) polarization for left boundary configurations
# CompuCellSetup.register_steppable(steppable=UniformGrowthSteppable(1,cd,tV,LamV,pP,pD,pC,growth_rate,t_relax)) #Implements cell growth for uniform cell proliferation 
# CompuCellSetup.register_steppable(steppable=MitosisUniformGrowthSteppable(1,cd,tV,LamV,pP,pD,pC,t_relax))       #Implements cell division for uniform cell proliferation
# CompuCellSetup.register_steppable(steppable=GrowthOnFrontSteppable(1,cd,tV,LamV,pP,pD,pC,growth_rate,t_relax))  #Implements cell growth for cell proliferation on front
# CompuCellSetup.register_steppable(steppable=MitosisOnFrontSteppable(1,cd,tV,LamV,pP,pD,pC,t_relax))             #Implements cell division for cell proliferation on front
CompuCellSetup.register_steppable(steppable=piff_generator(10000,Time,MCS0))                                    #Generates and stores .piff files
CompuCellSetup.register_steppable(steppable=adhesion_fix(10))                                                   #Ensures proximal and distal compartments are spatially segregated inside each cell
CompuCellSetup.run()
