"""
Code for Sub-cellular Potts Model for Planar Cell Polarity

Part of the paper
"Start Small: A Model for Tissue-wide Planar Cell Polarity without Morphogens"
Abhisha Thayambath and Julio M Belmonte

Code written by Abhisha Thayambath and Julio Monti Belmonte
Department of Physics, North Carolina State University

[Modified on 03/24/2025]

"""
from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
from math import *
from random import *
import time
from pathlib import Path
import os
import numpy as np
        

class InitialConditionsSteppable(SteppableBasePy):  
    """
    This Steppable initializes the simulation by creating an array of cells on the grid, with 
    specific target volumes for each cell-compartment type. The cells are composed of domains: PROX2, DIST2, CYTO and LAT. 
    The class sets up these domains by randomly distributing pixels from each CYTO cell into the different 
    domains based on defined pvolume fractions. The model starts with a DIST2 signal cell (first two columns) at the left boundary. 
    """
    def __init__(self,_frequency,_cd,_tV,_LamV,_pP,_pD,_pC):
        SteppableBasePy.__init__(self,_frequency)  
        # Initialize parameters for the steppable
        self.cd = _cd          # Cell diameter
        self.tV = _tV          # Target volume 
        self.LamV = _LamV      # Lambda volume 
        self.pP = _pP          # Volume fraction for PROX2 cells
        self.pD = _pD          # Volume fraction for DIST2 cells
        self.pC = _pC          # Volume fraction for LAT cells
        self.pL = 1 - _pP - _pD - _pC  # Remaining volume fraction for LAT cells (ensures total volume of a cell=self.tV)

    def start(self):
        # Create a DIST2 signal cell at the origin 
        signal = self.new_cell(self.DIST2)
        self.cellField[0, 0, 0] = signal
        
        # Create signal cells along the left boundary 
        for y in range(self.dim.y):
            x = 0
            self.cellField[x, y, 0] = signal
            x = 1
            self.cellField[x, y, 0] = signal
        
        # Set properties for the signal cell 
        signal.targetVolume = signal.volume
        signal.lambdaVolume = 1000
        
        # Drawing cells: Initialize CYTO cells in a grid-like pattern
        yy = self.cd / 2  # Starting y-coordinate
        for x0 in range(2, 9 * self.cd, self.cd):  # x-coordinates for placement of cells
            yy += self.cd / 2  # Adjust y-coordinate
            for y0 in range(0, self.dim.y, self.cd):  # y-coordinates for placement of cells
                x = x0
                y = (yy + y0) % self.dim.y  # Wrap around y-coordinate 
                cell = self.new_cell(self.CYTO)  # Create a new CYTO cell
                # Set properties for CYTO cells
                cell.targetVolume = self.tV
                cell.lambdaVolume = self.LamV
                # Place the CYTO cells in the simulation grid
                for x in range(self.cd):
                    x = x0 + x
                    for y in range(self.cd):
                        y = (yy + y0 + y) % self.dim.y  # Wrap around y-coordinate
                        self.cellField[x, y, 0] = cell
        
        # Setting cell properties for all CYTO cells
        for cell in self.cell_list_by_type(self.CYTO):
            cell.targetVolume = self.tV
            cell.lambdaVolume = self.LamV
        
        # get pixels associated with each cell
        self.get_all_pixels()
        
        for cell in self.cell_list_by_type(self.CYTO):
            # Assign pixels to their corresponding domains
            self.pixelate(cell)
    
    def pixelate(self, cell):
        # Create new cells for PROX2, DIST2, and LAT domains
        pCell = self.new_cell(self.PROX2)  # Create PROX2 cell
        dCell = self.new_cell(self.DIST2)  # Create DIST2 cell
        lCell = self.new_cell(self.LAT)    # Create LAT cell
        
        pList = []  # List to hold pixel locations

        # Add all pixels of all cells to pList
        for pixel in self.L[cell.id]:
            pList.append(pixel)
        
        shuffle(pList)  # Shuffle the list of pixels to randomize the pixel assignment
        
        # Assign pixels to different domains based on random probability
        for x, y, z in pList:
            r = random()  # Generate a random number
            if r < self.pP:
                self.cell_field[x, y, z] = pCell  # Assign to PROX2 domain
            elif r < self.pP + self.pD:
                self.cell_field[x, y, z] = dCell  # Assign to DIST2 domain
            elif r < self.pP + self.pD + self.pL:
                self.cell_field[x, y, z] = lCell  # Assign to LAT domain

        # Reassign the cluster IDs for the new domains (PROX2, DIST2, LAT) to that of the CYTO cell
        reassignIdFlag = self.reassign_cluster_id(pCell, cell.clusterId)
        reassignIdFlag = self.reassign_cluster_id(dCell, cell.clusterId)
        reassignIdFlag = self.reassign_cluster_id(lCell, cell.clusterId)
        
    def get_all_pixels(self):
        # Initialize a dictionary to store the pixels of each cell
        self.L = {}
        
        # Iterate over all cells and initialize a list for each cell
        for cell in self.cell_list:
            self.L[cell.id] = []
        
        # Iterate over all pixels in the grid
        for x, y, z in self.every_pixel():
            cell = self.cell_field[x, y, 0]  # Get the cell at the current pixel
            if cell:
                # Add the pixel coordinates to the list for the corresponding cell
                self.L[cell.id].append([x, y, 0])
 

#METRICS:  
class Volume_steppable(SteppableBasePy): #Sets volume constraint for cells when using Piff Initializer
    """
    This Steppable enforces volume constraints on different cell types in the simulation.

    """
    def __init__(self,_frequency,_tV,_LamV):
        SteppableBasePy.__init__(self,_frequency)
        self.tV=_tV; self.LamV=_LamV;
        
    def start(self):          
        self.set_default_volume()       
                   
    def set_default_volume(self):
        for cell in self.cell_list_by_type(self.CYTO):
            cell.targetVolume=cell.volume        # Maintain initial volume
            cell.lambdaVolume=self.LamV*2
        for cell in self.cell_list_by_type(self.PROX2):
            cell.targetVolume=cell.volume
            cell.lambdaVolume=self.LamV*2 
        for cell in self.cell_list_by_type(self.DIST2):
            # cell.targetVolume=cell.volume                 #for periodic boundary configuration
            # cell.lambdaVolume=self.LamV*2
            if cell.xCOM<2:                                     #for left boundary signal configuration
                cell.targetVolume=cell.volume
                cell.lambdaVolume=1000
            else:
                cell.targetVolume=cell.volume
                cell.lambdaVolume=self.LamV*2
        for cell in self.cell_list_by_type(self.LAT):
            cell.targetVolume=cell.volume
            cell.lambdaVolume=self.LamV*2    

            
                         
            
class CellPol(SteppableBasePy):
    """
    This Steppable computes the polarity of cells based on the relative positions
    of their proximal (PROX2) and distal (DIST2) components.

    Arguments:
          _frequency -- Frequency at which this Steppable runs
    """

    def __init__(self,_frequency,_cd):
        SteppableBasePy.__init__(self,_frequency)
        self.cd = _cd

    def start(self):
        for cell in self.cell_list_by_type(self.CYTO): 
            cell.dict["Ang-x"]=0; cell.dict["Ang-y"]=0  #Sets the x and y components of initial polarity vectors to zero    
    
    def step(self,mcs):
        for cell in self.cell_list_by_type(self.CYTO): 
            compList=self.get_cluster_cells(cell.clusterId)     # Get all cells in the cluster
            xP=0.0; yP=0.0; xD=0.0; yD=0.0; ok=0; p=0           # Initialize variables
            
            # Iterate through cluster cells to find the center of mass of proximal and distal compartments
            for cell2 in compList:
                if (cell2.type==self.PROX2): #proximal
                    xP=cell2.xCOM; yP=cell2.yCOM; ok+=1
                elif (cell2.type==self.DIST2): #distal
                    xD=cell2.xCOM; yD=cell2.yCOM; ok+=1
                    if (cell2.type==self.DIST2):
                        p+=1
            if (ok==2):      # If both proximal and distal components are found, compute polarity
                dx=xP-xD; dy=yP-yD  # Difference in center of mass of x and y-coordinates
                
                # Apply periodic boundary conditions
                if (abs(dy)>self.dim.y/2): dy=-(self.dim.y-abs(dy))*dy/abs(dy)
                if (abs(dx)>self.dim.x/2): dx=-(self.dim.x-abs(dx))*dx/abs(dx)
                
                d=sqrt(dx*dx+dy*dy) # Compute Euclidean distance between proximal and distal compartments
                
                if (d>0):
                    cell.dict["Ang-x"]=dx/d     # Normalize x-component of polarity
                    cell.dict["Ang-y"]=dy/d     # Normalize y-component of polarity
                    cell.dict["Pol"]=d/self.cd  # Normalize polarity magnitude
                else:
                    # If distance is zero, set polarity to zero
                    cell.dict["Ang-x"]=0; cell.dict["Ang-y"]=0; cell.dict["Pol"]=0
            else:
                # If proximal or distal component is missing, set polarity to zero
                cell.dict["Ang-x"]=0; cell.dict["Ang-y"]=0; cell.dict["Pol"]=0
            
            

class GraphGlobalPol(SteppableBasePy):
    """
    This Steppable calculates the global polarization (Phi) of cells over time and 
    stores the results in an output file. It also tracks the angular distribution of polarity vectors.
    
    Arguments:
          _frequency -- Frequency at which this Steppable runs
          _cd -- cell diameter used for normalizing polarity
          _t_piff -- Time offset for writing data when Piff_initializer is used

    """
    def __init__(self, _frequency, _cd, _t_piff, _growth_criteria):
        SteppableBasePy.__init__(self, _frequency)
        self.cd = _cd
        self.t_piff = _t_piff
        self.growth_criteria = _growth_criteria
        
    def start(self):
        # Initializes data structures for storing results.
        self.M = []     # Stores Monte Carlo steps (MCS)
        self.Phi_matrix = []    # Stores global polarization values
        self.angle = []         # Stores mean angles of polarization
        self.n_cell_matrix = [] # Stores the number of clusters/cells in the system
        self.inner_product_data = []
        
        # Initialize polarity values
        for cell in self.cell_list_by_type(self.CYTO): #Cytosol
            cell.dict["Ang-x"] = 0
            cell.dict["Ang-y"] = 0          
        
    def step(self, mcs): 
        if mcs >= 0:
            n = 0       # Counter for cells
            x = 0.0     # Sum of x-components of polarity
            y = 0.0     # Sum of y-components of polarity
            
            inner_pdt_data = 0
            count = 0

            # Compute average polarization vector
            for cell in self.cell_list_by_type(self.CYTO):                
                n += 1
                x += cell.dict['Ang-x']
                y += cell.dict['Ang-y']
                
                if self.growth_criteria == 0:
                    cell_id_list_for_growth = [34,35,36,37,153,157,165,169,177,181,189,193] #FOR GROWTH ON FRONT
                elif self.growth_criteria == 2:
                    cell_id_list_for_growth = [20, 153, 157]                        # FOR SINGLE CELL
                elif self.growth_criteria == 1:
                    cell_id_list_for_growth = [18,19,20,21,153,157,165,169,177,181,189,193] #FOR COLUMN OF CELLS INSIDE

                if cell.id in cell_id_list_for_growth:  #FOR COLUMN OF CELLS INSIDE               
                    inner_pdt_data += cell.dict["Ang-x"]
                    count += 1

            if count > 0:
                inner_pdt_data /= count

            self.inner_product_data.append(inner_pdt_data) 
            Phi = 1 # Default global polarization
            
            
            
            if n > 0:
                Phi = sqrt(x*x + y*y) / n       # Calculate average polarization by normalizing by cell count
                meanAng = pi + pi/2 - atan2(x, y) # Compute mean angle of polarization
               
                # Ensure angle is within [0, 360) degrees
                if degrees(meanAng) > 360:
                    angle = degrees(meanAng) - 360
                else:
                    angle = degrees(meanAng)

                # Store angle in the range [-180,180]
                if 0 <= angle < 180:
                    self.angle.append(angle)
                elif 180 <= angle <= 360:
                    self.angle.append(angle - 360)

                # Compute angular standard deviation
                std = 0.0
                nn = 0      
                for cell in self.cell_list_by_type(self.CYTO): #Cytosol
                    xx = cell.dict["Ang-x"]
                    yy = cell.dict["Ang-y"]
                    if xx and yy:
                        d = pi + pi/2 - atan2(xx, yy) - meanAng
                        if d > +pi:
                            d -= pi
                        if d < -pi:
                            d += pi
                        std += d*d
                        nn += 1

            # Count the number of clusters/cells in the simulation                    
            n_clusters = 0
            for compartments in self.clusters:        
                n_clusters += 1

            # Store the computed values
            self.n_cell_matrix.append(n_clusters)
            self.M.append(mcs + self.t_piff)
            self.Phi_matrix.append(Phi)     
            
            self.shared_steppable_vars['global_pol'] = Phi
            
            # Save results to an output file
            if self.output_dir is not None:
                 # create folder to store data
                output_path = Path(self.output_dir).joinpath("globalpol.txt")

                 # Ensure the output directory exists
                output_path.parent.mkdir(parents=True, exist_ok=True)
                try:
                    file_handle = open(output_path, 'a')
                    file_handle.write('{} {} {} {} {}\n'.format(mcs + self.t_piff, Phi, n_clusters, self.angle[-1], -1*self.inner_product_data[-1]))
                except IOError:
                    print ("Could not open file for writing.")
                    return
        


class piff_generator(SteppableBasePy):
    """
    Generates and stores .piff files
    """
    def __init__(self,_frequency, time, _t_piff):
        SteppableBasePy.__init__(self,_frequency)
        self.time=time; self.t_piff=_t_piff
        
    def step(self,mcs):
        self.get_all_pixels()               # Retrieve all cell pixels without using a pixel tracker
        self.save_piff(self.t_piff+mcs)                # Generate and save the PIF file
        
    def finish(self):
        self.get_all_pixels()               # Retrieve all cell pixels without using a pixel tracker
        self.save_piff(self.time+self.t_piff)               # Generate and save the PIF file
                    
                    
    def save_piff(self, mcs):
        out_folder = self.output_dir
        FileName = out_folder+"/PiffFile_"+str(mcs)+".piff"
        piffPath=Path(out_folder).joinpath(FileName)

        with open(piffPath, 'a') as fout:
            fout.write("Include Clusters \n")

            # Iterate through all cells stored in dictionary L
            for i in self.L:
                pixel_list = self.L[i]            # Get the list of pixels for the cell    
                cell = self.fetch_cell_by_id(i)   
                name = self.get_type_name_by_cell(cell)

                # Write pixel data in the PIF format
                for pixel in pixel_list:
                    x = pixel[0]; y = pixel[1]; z = pixel[2]          
                    fout.write("%d %d %s %d %d %d %d %d %d \n" % (cell.clusterId, cell.id, name, x, x, y, y, z, z))
    
    def get_all_pixels(self):
        self.L = {}   # Initialize a dictionary to store cell pixels
        for cell in self.cell_list:
            self.L[cell.id] = []

        # Iterate over all pixels in the simulation domain
        for x, y, z in self.every_pixel():
            cell = self.cell_field[x, y, 0]       # Get the cell occupying the pixel
            if cell:                             # If a cell exists at the given coordinates
                self.L[cell.id].append([x,y,0])     # Store the pixel coordinates
       


class Growth(SteppableBasePy):
    """ This class implements growth for cell proliferation on front.
    Arguments: 
            frequency -- Frequency at which this Steppable runs
                _cd -- cell diameter   
                _tV -- target volume of the cell
                _LamV -- lambda volume of the cell  
                _pP -- proportion of proximal domain
                _pD -- proportion of distal domain
                _pC -- proportion of cytoplasm domain
                _t_relax -- relaxation time for division
    """
    def __init__(self, _frequency, _cd,_tV,_LamV,_pP,_pD,_pC,_growth_rate,t_relax, _growth_criteria):
        SteppableBasePy.__init__(self, _frequency)
        self.cd=_cd; self.tV=_tV; self.LamV=_LamV; 
        self.pP=_pP; self.pD=_pD; self.pC=_pC; self.pL=1-_pP-_pD-_pC
        self.growth_rate = _growth_rate
        self.t_relax = t_relax
        self.growth_criteria = _growth_criteria
    
    
    def step(self, mcs):            
        if self.shared_steppable_vars['global_pol'] > 0.9:      # Growth starts only after relaxation time has passed
            growth_cluster_id_list = []                     # list to store cluster ids of cells to grow
            shrink_cluster_id_list = []                     # list to store cluster ids of cells to shrink
            cells_to_shrink = shrink_cluster_id_list  
            if self.growth_criteria == 0:
            #*********************************GROW ON FRONT*************************************************#
                for cell in self.cell_list_by_type(self.CYTO):      # goes over all clusters
                    if (self.cd<=cell.xCOM<=(self.dim.x - 2*self.cd)): # check if the cell is within the dimensions of the lattice       
                        # calculate the volume of the cluster
                        cluster_volume = 0       
                        for cell_cmpt in self.get_cluster_cells(cell.clusterId):   
                            cluster_volume += cell_cmpt.targetVolume  
                            # check if the cell is on the front (by checking if the cell/compartment has no neighbors)
                            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell_cmpt):
                                if (not neighbor):
                                    if not cell.dict.get("divided", False):
                                        growth_cluster_id_list.append(cell.clusterId)
                                        break
                        # check if the cell is not on the front and has a volume greater than target volume
                        if (cell.clusterId not in growth_cluster_id_list) and (cluster_volume>self.tV):
                            shrink_cluster_id_list.append(cell.clusterId)   # save the cluster id of the cell to shrink
                cells_to_grow = np.unique(growth_cluster_id_list)   # removes duplicate cluster ids
                cells_to_shrink = shrink_cluster_id_list
            #************************************************************************************************#                           
            
            
                      
            elif self.growth_criteria == 1:  
            #***********************************GROW ONE COLUMN OF CELLS INSIDE*************************************#
                cells_to_grow = [18, 19, 20, 21]
            #************************************************************************************************#
            elif self.growth_criteria == 2:
            #***********************************GROW ONE CELL INSIDE*************************************#
                cells_to_grow = [20]
            #************************************************************************************************#
            
            
            # filter out clusters where a CYTO cell has divided == 1
            cells_to_grow = [
                cluster_id
                for cluster_id in cells_to_grow
                if not any(
                    cell.type == self.CYTO and cell.dict.get("divided", 0) == 1
                    for cell in self.get_cluster_cells(int(cluster_id))
                )
            ]
     
            # if (mcs==self.t_relax):                                                                   #freeze the cells/clusters other than the ones growing
                # for cell in self.cell_list_by_type(self.CYTO):
                    # if cell.clusterId not in cells_to_grow:
                        # for cell_cmpt in self.get_cluster_cells(cell.clusterId):
                            # cell_cmpt.targetVolume=cell_cmpt.volume
                            # cell_cmpt.lambdaVolume=1000
            
            for cluster_to_grow in cells_to_grow:   # going over all cells to grow
                cluster_volume = 0          
                for cell in self.get_cluster_cells(int(cluster_to_grow)):
                    cluster_volume += cell.targetVolume     # estimate the current target volume of the cells
                    
                    
                if(cluster_volume<2*self.tV):       # if the current target volume is less than twice the target volume, increase the volume of the cell
                    # adding volume
                    cluster_volume+=self.growth_rate
                    S= 4*np.sqrt(cluster_volume)
                    
                    for cell in self.get_cluster_cells(int(cluster_to_grow)):                       # assign the target volume to each cell compartment   
                        if(cell.type==self.CYTO):
                            cell.targetVolume = (cluster_volume - 1.5*S)
                        elif (cell.type==self.PROX2):
                            cell.targetVolume = 0.45*S
                        elif (cell.type==self.DIST2):
                            cell.targetVolume = 0.45*S
                        elif (cell.type==self.LAT):
                            cell.targetVolume = 0.6*S  
                            
            for cluster_to_shrink in cells_to_shrink:   # going over all cells to shrink
                cluster_volume = 0          
                for cell in self.get_cluster_cells(int(cluster_to_shrink)):
                    cluster_volume += cell.targetVolume
                if(cluster_volume>self.tV):                    # if the current target volume is greater than the target volume, decrease the volume of the cell
                    # removing volume
                    cluster_volume-=self.growth_rate
                    S= 4*np.sqrt(cluster_volume)
                    
                    for cell in self.get_cluster_cells(int(cluster_to_shrink)):                        # assign the target volume to each cell compartment
                        if(cell.type==self.CYTO):
                            cell.targetVolume = (cluster_volume - 1.5*S)
                        elif (cell.type==self.PROX2):
                            cell.targetVolume = 0.45*S
                        elif (cell.type==self.DIST2):
                            cell.targetVolume = 0.45*S
                        elif (cell.type==self.LAT):
                            cell.targetVolume = 0.6*S 
                            


class Mitosis(MitosisSteppableClustersBase):
    """ This class implements mitosis for cell proliferation on front.
    Arguments: 
            frequency -- Frequency at which this Steppable runs
            _cd -- cell diameter 
            _tV -- target volume of the cell
            _LamV -- lambda volume of the cell
            _pP -- proportion of proximal domain
            _pD -- proportion of distal domain
            _pC -- proportion of cytoplasm domain
            _t_relax -- relaxation time for division
            
    """
    def __init__(self, frequency,_cd,_tV,_LamV,_pP,_pD,_pC,t_relax,_axis_to_divide):
        MitosisSteppableClustersBase.__init__(self, frequency)
        self.cd=_cd; self.tV=_tV; self.LamV=_LamV; 
        self.pP=_pP; self.pD=_pD; self.pC=_pC; self.pL=1-_pP-_pD-_pC
        self.t_relax = t_relax
        self.axis_to_divide = _axis_to_divide
    
    def start(self):
        self.time_of_division = []
    
    def step(self, mcs):
        
        cells_to_divide = []
        if(mcs>=0):      # Growth starts only after relaxation time has passed
            mitosis_cluster_id_list = []
            for cell in self.cell_list_by_type(self.CYTO):
                cluster_cell_list = self.get_cluster_cells(cell.clusterId)
                cluster_volume = 0
                for cell_cmpt in cluster_cell_list:
                    cluster_volume += cell_cmpt.volume
                if(cluster_volume>=2*self.tV):# :       # if the current target volume is greater than twice the target volume, divide the cell
                    mitosis_cluster_id_list.append(cell.clusterId) 
            for cluster_to_divide in mitosis_cluster_id_list:           # going over all cells to divide
                print("Time of division", mcs, cluster_to_divide,"dividing - yes")    
                # self.divide_cluster_random_orientation(cluster_to_divide)   #to divide the cell along random axis
                if self.axis_to_divide == 0:
                    self.divide_cluster_along_minor_axis(cluster_to_divide)         # to divide the cell along minor axis
                elif self.axis_to_divide == 1:
                    self.divide_cluster_along_major_axis(cluster_to_divide)         # to divide the cell along minor axis
                #self.divide_cluster_along_major_axis(cluster_to_divide)         # to divide the cell along major axis
                self.check()        # to check if the cell is divided and inherited compartments correctly
                self.check_for_ratios() # to check if the cell is divided and inherited compartments correctly
                self.time_of_division.append(mcs)
                

                # Save results to an output file
                if self.output_dir is not None:
                     # create folder to store data
                    output_path = Path(self.output_dir).joinpath("time_of_division.txt")

                     # Ensure the output directory exists
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    try:
                        file_handle = open(output_path, 'a')
                        file_handle.write('{}\n'.format(mcs))
                    except IOError:
                        print ("Could not open file for writing.")
                        return
    
    def check(self):    # to fix the bug that doesn't assign the correct target volume to the newly divided cell
        #self.cells_to_delete = []
        for cell in self.cell_list_by_type(self.CYTO):  
            cluster_cell_list = self.get_cluster_cells(cell.clusterId)
            cluster_volume = 0
            for cell_cmpt in cluster_cell_list:
                cluster_volume += cell_cmpt.targetVolume
            if cluster_volume > (2*self.tV + 0.5):          # check if the target volume of the cell is greater than twice the target volume even after division
                print("#########################CELL IS WRONG##########################")
                print("cluster volume",cluster_volume)
                for cell_cmpt in cluster_cell_list:
                    print(cell_cmpt.type, " volume " ,cell_cmpt.volume," target volume ", cell_cmpt.targetVolume," lambda volume ", cell.lambdaVolume)
                    # assign the target volume to each cell compartment of the newly divided cell which did not update its target volumes
                    if(cell_cmpt.type==self.PROX2):
                        cell_cmpt.targetVolume=(self.tV)*self.pP
                    elif(cell_cmpt.type==self.DIST2):
                        cell_cmpt.targetVolume=(self.tV)*self.pD
                    elif(cell_cmpt.type==self.LAT):
                        if(cell_cmpt.volume==0):
                            cell_cmpt.targetVolume = 0
                            #self.cells_to_delete.append(cell_cmpt.id)
                        else:
                            cell_cmpt.targetVolume=(self.tV)*self.pL
                    elif(cell_cmpt.type==self.CYTO):
                        cell_cmpt.targetVolume=(self.tV)*(1-self.pL-self.pP-self.pD)
                    print("%%%%%%new volume%%%%%%")
                    print(cell_cmpt.type, " volume " ,cell_cmpt.volume," target volume ", cell_cmpt.targetVolume," lambda volume ", cell.lambdaVolume)
        #for cell in self.cells_to_delete:
            #self.delete_cell(cell)
            
    def check_for_ratios(self):     # to fix the bug that doesn't assign the correct target volume to the newly divided cell
        for cell in self.cell_list_by_type(self.CYTO):
            cluster_cell_list = self.get_cluster_cells(cell.clusterId)
            cluster_volume = 0
            for cell_cmpt in cluster_cell_list:
                cluster_volume += cell_cmpt.targetVolume
                if(cell_cmpt.type==self.CYTO):
                    cyto_vol = cell_cmpt.volume
                    cyto_target_vol = cell_cmpt.targetVolume
            if cluster_volume > (2*self.tV + 0.5):
                print("#########################CELL IS WRONG##########################")
                print("cluster volume",cluster_volume)
                for cell_cmpt in cluster_cell_list:
                    print(cell_cmpt.type, " volume ratio" ,cell_cmpt.volume/cyto_vol," target volume ratio", cell_cmpt.targetVolume/cyto_target_vol)
            
    
    def update_attributes(self):      
        self.compartment_list_parent = self.get_cluster_cells(self.parent_cell.clusterId)       # get the list of compartments of the parent cell
        self.compartment_list_child = self.get_cluster_cells(self.child_cell.clusterId)         # get the list of compartments of the child cell
        print("printiing ids of parent and daughter cells", self.parent_cell.id,self.child_cell.id)
        self.pixelate(self.compartment_list_parent, self.parent_cell)   # pixelate the parent cell (uniformly mix all the inherited domains)
        self.pixelate(self.compartment_list_child, self.child_cell)     # pixelate the child cell (uniformly mix all the inherited domains)

        
    def pixelate(self, compartment_list, cell):
        pCell = self.new_cell(self.PROX2)            #create proximal domain
        dCell = self.new_cell(self.DIST2)            #create distal domain
        lCell = self.new_cell(self.LAT)             #create lateral domain
        cytoCell = self.new_cell(self.CYTO)             #create lateral domain
        pList = []
        # Iterate through all pixels in the simulation domain
        for x, y, z in self.every_pixel():
            new_cell = self.cell_field[x, y, 0]       
            if new_cell:
                for cell_ in compartment_list:
                    if (new_cell.id==cell_.id):     #   check if the cell is one of the compartments of the newly divided cell
                        pList.append([x,y,0])

        shuffle(pList)      # shuffle the list of pixels of the newly divided cell
        cluster_volume = len(pList)
        p_frac = self.pP
        d_frac = self.pD
        l_frac = self.pL
  
        for x,y,z in pList:
            # randomly assign the cell type to each pixel of the newly divided cell
            r=random()
            if (r<p_frac):
                self.cell_field[x,y,z] = pCell
            elif (r<p_frac+d_frac):
                self.cell_field[x,y,z] = dCell
            elif (r<p_frac+d_frac+l_frac):
                self.cell_field[x,y,z] = lCell
            elif(r<1):
                self.cell_field[x,y,z] = cytoCell
        
    #CYTO subcell        
        cytoCell.targetVolume=(self.tV)*(1-p_frac-d_frac-l_frac)        # assign the target volume to each cell compartment
        cytoCell.lambdaVolume=self.LamV*2        
        reassignIdFlag=self.reassign_cluster_id(cytoCell,cell.clusterId) # reassign the cluster id to the newly divided cell
        cytoCell.dict["divided"] = 1
                
    #PROX subcell
        pCell.targetVolume=(self.tV)*p_frac
        pCell.lambdaVolume=self.LamV*2    
        reassignIdFlag=self.reassign_cluster_id(pCell,cell.clusterId)
    #DIST subcell
        dCell.targetVolume=(self.tV)*d_frac
        dCell.lambdaVolume=self.LamV*2
        reassignIdFlag=self.reassign_cluster_id(dCell,cell.clusterId)
    #LAT subcell
        lCell.targetVolume=(self.tV)*l_frac
        lCell.lambdaVolume=self.LamV*2
        reassignIdFlag=self.reassign_cluster_id(lCell,cell.clusterId)  
        
        print("printing ids of child cell", cytoCell.id)

class adhesion_fix(SteppableBasePy):
    """
    This class implements a fix for the adhesion of proximal-distal compartments (to prevent these from being nearest neighbors) inside each cell
    """
    def __init__(self,_frequency):
        SteppableBasePy.__init__(self,_frequency)
    
    def step(self,mcs):
        # Iterate through all cells in the simulation
        for cell in self.cell_list_by_type(self.CYTO):
            cluster_cell_list = self.get_cluster_cells(cell.clusterId)
            l = 0
            p = 0
            # Iterate through all compartments in the cluster
            for cell_cmpt in cluster_cell_list:
                if(cell_cmpt.type==self.PROX2):
                    p = cell_cmpt               # get the proximal compartment
                elif (cell_cmpt.type==self.LAT):
                    l = cell_cmpt               # get the lateral compartment
            if (l and p):   # If both proximal and lateral compartments are found
                # Get the list of boundary pixels for the proximal compartment
                pixel_list = self.get_cell_boundary_pixel_list(p)

                if pixel_list is not None:         
                    for boundary_pixel_tracker_data in pixel_list:  # Iterate through all boundary pixels
                        # Get the coordinates of the boundary pixels
                        x = boundary_pixel_tracker_data.pixel.x 
                        y = boundary_pixel_tracker_data.pixel.y
                        z = boundary_pixel_tracker_data.pixel.z

                        # Get the nearest neighboring cells of the boundary pixels
                        neighbors = [(x + dx, y + dy) for dx in [-1, 0, 1] for dy in [-1, 0, 1]  if (dx, dy) != (0, 0)]
                        
                        # Check if any neighboring cell is of type DIST2 and belongs to the same cluster
                        for nx, ny in neighbors:
                            neighbor_cell = self.cell_field[nx, ny, 0]
                            if neighbor_cell:
                                if (neighbor_cell.type == self.DIST2 and neighbor_cell.clusterId == cell.clusterId):
                                    # If found, set the boundary pixel to the lateral compartment
                                    self.cell_field[x, y, z] = l
                                    break
                                    
