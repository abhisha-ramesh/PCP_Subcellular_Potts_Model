# Code written by Abhisha Thayambath and Julio Monti Belmonte
# [Department of Physics, North Carolina State University]
# [Modified on 03/24/2025]
global cd, Lx, Ly, T, nOrder, Time
global LamV, tV
global pP, pD, pC
global E_pd, E_dd, E_pp, E_ll, E_lp, E_ld, E_ml, J_pd, E_pd_mut, J_pd_mut, E_dmut

#PARAMETERS:
cd=12           #Typical cell diameter
#Lx=8*cd        #Size of lattice (Lx X Ly) for periodic boundary conditions
Lx=(9*cd)+cd   #Size of lattice (Lx X Ly) for left boundary signal
Ly=8*cd        #Size of lattice (Lx X Ly)
#
#POTTS PARAMETERS:
T=15            #Temperature
nOrder=4        #Distance of interaction
Time=100000    #Total number of MCS 
#
#VOLUME/SURFACE PARAMETERS:
LamV=12         #Lambda Volume
tV=cd*cd        #Target Volume
#
#SUBCELLULAR PARAMETERS:
pP=0.15        #Proximal domain
pD=0.15        #Distal domain
pC=0.50        #Cytoplamsic domain
#

#CONTACT ENERGIES:
E_pd = 7                #prox-dist
E_dd = 12              #dist-dist
E_lp = 12               #lat-prox
E_ll = 14               #lat-lat
E_ld = E_lp             #lat-dist
E_ml = E_ll             #medium-lat
E_pp = E_dd             #prox-prox


E_pd_mut = 7
E_dmut = 12
#
#INTERNAL CONTACT ENERGIES         #Internal contact energy between proximal and distal compartments
J_pd = 30 
J_pd_mut = 30
#

def configure_simulation():

    from cc3d.core.XMLUtils import ElementCC3D
    
    CompuCell3DElmnt=ElementCC3D("CompuCell3D",{"Revision":"0","Version":"4.3.1"})
    
    MetadataElmnt=CompuCell3DElmnt.ElementCC3D("Metadata")
    
    MetadataElmnt.ElementCC3D("NumberOfProcessors",{},"1")
    MetadataElmnt.ElementCC3D("DebugOutputFrequency",{},"0")
    
    PottsElmnt=CompuCell3DElmnt.ElementCC3D("Potts")
    
    PottsElmnt.ElementCC3D("Dimensions",{"x":Lx,"y":Ly,"z":"1"})
    PottsElmnt.ElementCC3D("Steps",{},Time)
    PottsElmnt.ElementCC3D("Temperature",{},int(T))
    PottsElmnt.ElementCC3D("NeighborOrder",{},nOrder)
    PottsElmnt.ElementCC3D("Boundary_x",{},"Periodic")
    PottsElmnt.ElementCC3D("Boundary_y",{},"Periodic")
    
    
    
    PluginElmnt=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CellType"})
    
    # Listing all cell types in the simulation
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"0","TypeName":"Medium"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"1","TypeName":"cyto"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"2","TypeName":"prox"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"3","TypeName":"dist"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"4","TypeName":"lat"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"5","TypeName":"mut"})
    
    CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"VolumeLocalFlex"})
    
    PluginElmnt_1=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CenterOfMass"})
    
    PluginElmnt_2=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"NeighborTracker"})
    PluginElmnt_7=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"PixelTracker"})
    PluginElmnt_4=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Contact"})
    
    
    # Specification of adhesion energies
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"Medium"},"0.0")
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"cyto"},"30")
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"prox"},"10.0")
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"dist"},"10.0")
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"lat"},E_ml)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"Medium","Type2":"mut"},"10")
    
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"cyto","Type2":"cyto"},"30")
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"cyto","Type2":"prox"},"30")
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"cyto","Type2":"dist"},"30")
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"cyto","Type2":"lat"},"30")    
    
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"prox","Type2":"prox"},E_pp)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"prox","Type2":"dist"},E_pd)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"prox","Type2":"lat"},E_lp)

    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"dist","Type2":"dist"},E_dd)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"dist","Type2":"lat"},E_ld)
    
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"lat","Type2":"lat"},E_ll)
    
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"mut","Type2":"mut"},E_dd)    #because mutant cell type is distal(Vang mutant)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"mut","Type2":"dist"},E_dmut)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"mut","Type2":"prox"},E_pd_mut)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"mut","Type2":"lat"},E_ld)
    PluginElmnt_4.ElementCC3D("Energy",{"Type1":"cyto","Type2":"mut"},"30")
    
    PluginElmnt_4.ElementCC3D("NeighborOrder",{},nOrder)
    
    
    PluginElmnt_5=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"ContactInternal"})
    # Specification of internal adhesion energies
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"cyto","Type2":"prox"},"1.0")
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"cyto","Type2":"mut"},"1.0")
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"cyto","Type2":"dist"},"1.0")
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"cyto","Type2":"lat"},"1.0")
    
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"prox","Type2":"dist"},J_pd)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"prox","Type2":"lat"},"1.0")
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"mut","Type2":"prox"},J_pd_mut)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"mut","Type2":"lat"},"1.0")

    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"dist","Type2":"lat"},"1.0")
    PluginElmnt_5.ElementCC3D("NeighborOrder",{},nOrder)
    
    
    
    CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)
            
from cc3d import CompuCellSetup
        

configure_simulation()            

            

from PCP_simpleSteppables import *

CompuCellSetup.register_steppable(steppable=InitCondY_bXl(10,cd,tV,LamV,pP,pD,pC))  
CompuCellSetup.register_steppable(steppable=MutantPhenotype(10,cd,tV,LamV,pP,pD,pC))
CompuCellSetup.register_steppable(steppable=CellPol(10,cd))
CompuCellSetup.register_steppable(steppable=GraphGlobalPol(10,cd))
CompuCellSetup.register_steppable(steppable=OrientationField(10))
CompuCellSetup.register_steppable(steppable=VectorField(1000))
CompuCellSetup.register_steppable(steppable=piff_generator(1000,Time,0))
CompuCellSetup.run()
