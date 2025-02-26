from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
from math import *
from random import *
import time
from pathlib import Path
import os
import numpy as np
        
class periodic_boundary(SteppableBasePy):  
    """
    This Steppable defines the initial conditions for the periodic boundary configuration. Seeds of cells of type 
    "TEST" are randomly placed on the grid. A relaxation period is included to increase the volumes of the TEST cells to match 
    the target volume of a cell (self.tV), and the cells are later compartmentalised into different domains (PROX2, DIST2, CYTO and LAT) 
    after the relaxation time is passed.
    """
    
    def __init__(self,_frequency,_cd,_tV,_LamV,_pP,_pD,_pC,_nx,_ny,_t_relax):
        """
        Initialize the periodic boundary coconfiguration with necessary parameters.

        Arguments:
          _frequency -- Frequency at which this Steppable is executed
          _cd -- Cell diameter
          _tV -- Target volume for cells
          _LamV -- Lambda volume for cells
          _pP -- Volume fraction for PROX2 cells
          _pD -- Volume fraction for DIST2 cells
          _pC -- Volume fraction for LAT cells
          _nx -- Number of cells in the x-direction
          _ny -- Number of cells in the y-direction
          _t_relax -- Relaxation time for adjusting cell volume and relaxing cell boundaries
        """
        SteppableBasePy.__init__(self,_frequency)
        self.cd=_cd  # Cell diameter
        self.tV=_tV  # Target volume for cells
        self.LamV=_LamV  # Lambda volume for cells
        self.pP=_pP  # Volume fraction for PROX2 cells
        self.pD=_pD  # Volume fraction for DIST2 cells
        self.pC=_pC  # Volume fraction for LAT cells
        self.pL=1-_pP-_pD-_pC  # Remaining volume fraction for LAT cells
        self.nx = _nx  # Number of cells in x-direction
        self.ny = _ny  # Number of cells in y-direction
        self.t_relax = _t_relax  # Relaxation time for adjusting cell volume and relaxing cell boundaries
    
    def left_right_n_cells(self, n_order, x_f, x_init, y_init):
        """
        Helper function to determine the number of left and right cells based on the order and initial positions.

        Arguments:
          n_order -- Order of the hexagon
          x_f -- Number of cells in the x-direction
          x_init -- Initial x-coordinate
          y_init -- Initial y-coordinate

        Returns:
          left_cells -- Number of left-side cells
          right_cells -- Number of right-side cells
        """
        if (y_init % 2 == 0):
            if (n_order % 2 != 0):    
                left_cells = (x_f // 2)
                right_cells = x_f - left_cells
            else:
                right_cells = x_f // 2
                left_cells = x_f - right_cells
        else:
            if (n_order % 2 == 0):    
                left_cells = (x_f // 2)
                right_cells = x_f - left_cells
            else:
                right_cells = x_f // 2
                left_cells = x_f - right_cells
                
        return left_cells, right_cells
    
    def draw(self, x_initial_list, y_initial_list, n_order):
        """
        Function to draw and place cells on the grid randomly
        Arguments:
          x_initial_list -- List of initial x-positions for the cells
          y_initial_list -- List of initial y-positions for the cells
          n_order -- order of the hexagon 
        """
        for y_initial in range(len(y_initial_list)):
            y_f = (n_order*2)+1  # Adjusted y-length
            x_f = n_order  # Number of cells in x-direction
            x_init = x_initial_list[y_initial]  # Starting x-coordinate
            cell = self.new_cell(self.TEST)  # Create a new TEST cell

            midpoint = ((x_init + ((n_order+1) / 2)))  # Calculate midpoint for alignment

            for y in range(y_initial_list[y_initial], y_initial_list[y_initial] + y_f):
                n_y = y_initial_list[y_initial] + n_order
                if (y <= n_y):
                    x_f += 1
                elif (y > n_y):
                    x_f -= 1
                
                left_cells, right_cells = self.left_right_n_cells(n_order, x_f, x_init, y_initial_list[y_initial])
                
                for x in range(int(midpoint - left_cells), int(midpoint + right_cells)):
                    if (0 <= x <= self.dim.x):  # Ensure cell is within bounds
                        self.cellField[x, y, 0] = cell  # Place the cell on the grid

            # Set target volume and lambda volume for the cell
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 10
    
    def start(self):
        """
        Initialize the grid by generating initial positions for cells and drawing them on the grid.

        This method starts by creating random initial positions for the cells and then draws the cells 
        using the `draw` method.
        """
        n_order = 2  # order for hexagon
        y_initial_list = [randint(0, self.dim.y) for _ in range(self.nx*self.ny)]  # Random y-coordinates
        x_initial_list = [randint(0, self.dim.x) for _ in range(self.nx*self.ny)]  # Random x-coordinates

        self.draw(x_initial_list, y_initial_list, n_order)  # Draw cells on the grid
    
    def step(self, mcs):
        """
        During each time step, cell volumes are increased until it reached self.tV, and cells are compartmentalised 
        based on relaxation time.

        """
        for cell in self.cell_list_by_type(self.TEST):
            if (mcs < self.t_relax):  # Relaxation phase: increase cell volume 
                if cell.volume < self.tV:
                    cell.targetVolume += 1
                    cell.lambdaVolume = 10
        if (mcs == self.t_relax):  # Once relaxation is complete (cells relax boundaries), start pixelation
            self.get_all_pixels()
            
            for cell in self.cell_list_by_type(self.TEST):
                cell.lambdaVolume = 10
                cell.type = self.CYTO  # Change cell type to CYTO
                self.pixelate(cell)  # Assign pixels to domains
                
    def pixelate(self, cell):
        """
        Assign pixels from the CYTO cell to different domains (PROX2, DIST2, LAT) based on defined volume fractions.

        Arguments:
          cell -- The CYTO cell to be pixelated
        """
        pCell = self.new_cell(self.PROX2)  # Create PROX2 cell
        dCell = self.new_cell(self.DIST2)  # Create DIST2 cell
        lCell = self.new_cell(self.LAT)    # Create LAT cell
        
        pList = []  # List of pixels for the cell

        # Add all pixels of the cell to pList
        for pixel in self.L[cell.id]:
            pList.append(pixel)
        shuffle(pList)  # Shuffle the pixel list to randomize assignment
        
        # Assign pixels to different domains based on random probability
        for x, y, z in pList:
            r = random()  # Generate random number for domain assignment
            if (r < self.pP):
                self.cell_field[x, y, z] = pCell  # Assign to PROX2 domain
            elif (r < self.pP + self.pD):
                self.cell_field[x, y, z] = dCell  # Assign to DIST2 domain
            elif (r < self.pP + self.pD + self.pL):
                self.cell_field[x, y, z] = lCell  # Assign to LAT domain

        # Set target volume and lambda volume for the new domains
        cell.targetVolume = self.tV * self.pC
        cell.lambdaVolume = self.LamV * 2

        pCell.targetVolume = self.tV * self.pP
        pCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(pCell, cell.clusterId)  # Reassign cluster ID for PROX2

        dCell.targetVolume = self.tV * self.pD
        dCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(dCell, cell.clusterId)  # Reassign cluster ID for DIST2

        lCell.targetVolume = self.tV * self.pL
        lCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(lCell, cell.clusterId)  # Reassign cluster ID for LAT
        
        
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

class left_boundary(SteppableBasePy):  # Left boundary signal, y-periodic
    """
    This Steppable defines the initial conditions for the left boundary signal configuration. Seeds of cells of type 
    "TEST" are randomly placed on the grid. A relaxation period is included to increase the volumes of the TEST cells to match 
    the target volume of a cell (self.tV), and the cells are later compartmentalised into different domains (PROX2, DIST2, CYTO and LAT) 
    after the relaxation time is passed. At relaxation time, the first two columns of pixels are made type DIST which acts as the left boundary orienting signal. 
    """

    def __init__(self, _frequency, _cd, _tV, _LamV, _pP, _pD, _pC, _nx, _ny, _t_relax):
        """
        Initialize the left boundary orienting signal configuration with necessary parameters.

        Arguments:
          _frequency -- Frequency at which this Steppable is executed
          _cd -- Cell diameter
          _tV -- Target volume for cells
          _LamV -- Lambda volume for cells
          _pP -- Volume fraction for PROX2 cells
          _pD -- Volume fraction for DIST2 cells
          _pC -- Volume fraction for LAT cells
          _nx -- Number of cells in the x-direction
          _ny -- Number of cells in the y-direction
          _t_relax -- Relaxation time for adjusting cell volume and relaxing cell boundaries
        """
        SteppableBasePy.__init__(self, _frequency)
        self.cd = _cd  # Cell diameter
        self.tV = _tV  # Target volume for cells
        self.LamV = _LamV  # Lambda volume for cells
        self.pP = _pP  # Volume fraction for PROX2 cells
        self.pD = _pD  # Volume fraction for DIST2 cells
        self.pC = _pC  # Volume fraction for LAT cells
        self.pL = 1 - _pP - _pD - _pC  # Remaining volume fraction for LAT cells
        self.nx = _nx  # Number of cells in x-direction
        self.ny = _ny  # Number of cells in y-direction
        self.t_relax = _t_relax  # Relaxation time for cell volume adjustment and relaxing cell boundaries

    def left_right_n_cells(self, n_order, x_f, x_init, y_init):
        """
        Helper function to determine the number of left and right cells based on the order and initial positions.

        Arguments:
          n_order -- order of the hexagon
          x_f -- Number of cells in the x-direction
          x_init -- Initial x-coordinate
          y_init -- Initial y-coordinate

        Returns:
          left_cells -- Number of left-side cells
          right_cells -- Number of right-side cells
        """
        if y_init % 2 == 0:
            if n_order % 2 != 0:
                left_cells = x_f // 2
                right_cells = x_f - left_cells
            else:
                right_cells = x_f // 2
                left_cells = x_f - right_cells
        else:
            if n_order % 2 == 0:
                left_cells = x_f // 2
                right_cells = x_f - left_cells
            else:
                right_cells = x_f // 2
                left_cells = x_f - right_cells

        return left_cells, right_cells

    def draw(self, x_initial_list, y_initial_list, n_order):
        """
        Function to draw and place cells on the grid randomly.

        Arguments:
          x_initial_list -- List of initial x-positions for the cells
          y_initial_list -- List of initial y-positions for the cells
          n_order -- order of the hexagon
        """
        for y_initial in range(len(y_initial_list)):
            y_f = (n_order * 2) + 1  # Adjusted y-length
            x_f = n_order  # Number of cells in x-direction
            x_init = x_initial_list[y_initial]  # Starting x-coordinate
            cell = self.new_cell(self.TEST)  # Create a new TEST cell

            midpoint = (x_init + ((n_order + 1) / 2))  # Calculate midpoint for alignment

            for y in range(y_initial_list[y_initial], y_initial_list[y_initial] + y_f):
                n_y = y_initial_list[y_initial] + n_order
                if y <= n_y:
                    x_f += 1
                else:
                    x_f -= 1

                left_cells, right_cells = self.left_right_n_cells(n_order, x_f, x_init, y_initial_list[y_initial])

                for x in range(int(midpoint - left_cells), int(midpoint + right_cells)):
                    if 0 <= x <= self.dim.x:  # Ensure cell is within bounds
                        self.cellField[x, y, 0] = cell  # Place the cell on the grid

            # Set target volume and lambda volume for the cell
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 10

    def start(self):
        """
        Initialize the grid by generating random initial positions for cells and drawing them.

        This method starts by creating random initial positions for the cells within the leftmost 
        third of the grid and then calls the `draw` function.
        """
        n_order = 2  # Number of cells to place per row
        y_initial_list = [randint(0, self.dim.y) for _ in range(self.nx * self.ny)]
        x_initial_list = [randint(2, int(self.dim.x / 3)) for _ in range(self.nx * self.ny)]        #we place the seeds only in 1/3rd of the lattice along x-direction since there should be open space on the front

        self.draw(x_initial_list, y_initial_list, n_order)  # Draw cells on the grid

    def step(self, mcs):
        """
        During each time step, cell volumes are adjusted, and a boundary signal 
        is introduced at the leftmost column.

        """
        for cell in self.cell_list_by_type(self.TEST):
            if mcs < self.t_relax:  # Relaxation phase: adjust cell volume and cells relax boundaries
                if cell.volume < self.tV:
                    cell.targetVolume += 1
                    cell.lambdaVolume = 10
        if mcs == self.t_relax:
            signal = self.new_cell(self.DIST2)  # Create a signal cell of type DIST2
            self.cellField[0, 0, 0] = signal  # Place signal at left boundary
            for y in range(self.dim.y):
                self.cellField[0, y, 0] = signal
                self.cellField[1, y, 0] = signal
            signal.targetVolume = signal.volume
            signal.lambdaVolume = 1000
            self.get_all_pixels()

            for cell in self.cell_list_by_type(self.TEST):
                cell.lambdaVolume = 10
                cell.type = self.CYTO   # Change cell type to CYTO
                self.pixelate(cell)     # Assign pixels to domains

    def pixelate(self, cell):
        """
        Assign pixels from the CYTO cell to different domains (PROX2, DIST2, LAT) based on defined volume fractions.

        Arguments:
          cell -- The CYTO cell to be pixelated
        """
        pCell = self.new_cell(self.PROX2)       # Create PROX2 cell
        dCell = self.new_cell(self.DIST2)       # Create DIST2 cell
        lCell = self.new_cell(self.LAT)         # Create LAT cell

        pList = []
        # Add all pixels of the cell to pList
        for pixel in self.L[cell.id]:
            pList.append(pixel)
        shuffle(pList)  # Shuffle the pixel list to randomize assignment

        # Assign pixels to different domains based on random probability
        for x, y, z in pList:
            r = random()    # Generate random number for domain assignment
            if r < self.pP:
                self.cell_field[x, y, z] = pCell    # Assign to PROX2 domain
            elif r < self.pP + self.pD:
                self.cell_field[x, y, z] = dCell    # Assign to DIST2 domain
            elif (r < self.pP + self.pD + self.pL):
                self.cell_field[x, y, z] = lCell    # Assign to LAT domain

        # Set target volume and lambda volume for the new domains
        cell.targetVolume = self.tV * self.pC
        cell.lambdaVolume = self.LamV * 2

        pCell.targetVolume = self.tV * self.pP
        pCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(pCell, cell.clusterId)  # Reassign cluster ID for PROX2

        dCell.targetVolume = self.tV * self.pD
        dCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(dCell, cell.clusterId)  # Reassign cluster ID for DIST2

        lCell.targetVolume = self.tV * self.pL
        lCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(lCell, cell.clusterId)  # Reassign cluster ID for LAT
        
        
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
        for x0 in range(2, 4 * self.cd, self.cd):  # x-coordinates for placement of cells
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
            cell.targetVolume=cell.volume                 #for periodic boundary configuration
            cell.lambdaVolume=self.LamV*2
            # if cell.xCOM<2:                                     #for left boundary signal configuration
                # cell.targetVolume=cell.volume
                # cell.lambdaVolume=1000
            # else:
                # cell.targetVolume=cell.volume
                # cell.lambdaVolume=self.LamV*2
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
    def __init__(self, _frequency, _cd, _t_piff):
        SteppableBasePy.__init__(self, _frequency)
        self.cd = _cd
        self.t_piff = _t_piff
        
    def start(self):
        # Initializes data structures for storing results.
        self.M = []     # Stores Monte Carlo steps (MCS)
        self.Phi_matrix = []    # Stores global polarization values
        self.angle = []         # Stores mean angles of polarization
        self.n_cell_matrix = [] # Stores the number of clusters/cells in the system
        
        # Initialize polarity values
        for cell in self.cell_list_by_type(self.CYTO): #Cytosol
            cell.dict["Ang-x"] = 0
            cell.dict["Ang-y"] = 0          
        
    def step(self, mcs): 
        if mcs >= 10000:
            n = 0       # Counter for cells
            x = 0.0     # Sum of x-components of polarity
            y = 0.0     # Sum of y-components of polarity
            
            # Compute average polarization vector
            for cell in self.cell_list_by_type(self.CYTO):                
                n += 1
                x += cell.dict['Ang-x']
                y += cell.dict['Ang-y']
            
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
            
            # Save results to an output file
            if self.output_dir is not None:
                 # create folder to store data
                output_path = Path(self.output_dir).joinpath("globalpol.txt")

                 # Ensure the output directory exists
                output_path.parent.mkdir(parents=True, exist_ok=True)
                try:
                    file_handle = open(output_path, 'a')
                    file_handle.write('{} {} {} {}\n'.format(mcs + self.t_piff, Phi, n_clusters, self.angle[-1]))
                except IOError:
                    print ("Could not open file for writing.")
                    return
        

class OrientationField(SteppableBasePy):
    def __init__(self,_frequency=1):
        SteppableBasePy.__init__(self,_frequency)
        # Create vector fields at the cell level for tracking orientation
        self.create_vector_field_cell_level_py("Orientation")
        
    def step(self,mcs):
        if (mcs>=10000):
            # Clear the Orientation fields at the beginning of each step
            field = self.field.Orientation.clear()
            
            for cell in self.cell_list_by_type(self.CYTO):  # Iterate over all CYTO cells
                if (cell.dict["Pol"]!=0):       # Only update cells with nonzero polarity
                    field = self.field.Orientation  # Access the Orientation field
                    # Set the orientation vector using cell polarity angles
                    field[cell] = [-cell.dict["Ang-x"], -cell.dict["Ang-y"], 0]
                    

class VectorField(SteppableBasePy):
    def __init__(self,_frequency=1):
        SteppableBasePy.__init__(self,_frequency)
        
    def step(self,mcs):
        # Lists to store vector field data
        self.X = []; self.Y = []; self.X_com = []; self.Y_com = []; 

        for cell in self.cell_list_by_type(self.CYTO):
            if (cell.dict["Pol"]!=0):   # Only process cells with nonzero polarity
                x = cell.dict['Ang-x']; y = cell.dict['Ang-y']  # Retrieve x and y-components of polarity
                
                # Store the componenets of polarity vectors along with cell center of masses
                self.Y.append(-y)
                self.X.append(-x)
                self.Y_com.append(cell.yCOM)
                self.X_com.append(cell.xCOM)
                
        # Open a file in the simulation output folder to store vector field data
        output_file, path__ = self.open_file_in_simulation_output_folder("vector_field_data_"+str(mcs)+".txt", mode='w')
        for i in range(len(self.X)):
            output_file.write('{} {} {} {}\n'.format(self.X_com[i],self.Y_com[i],self.X[i],self.Y[i]))                       
                
        
class GraphPol_x_Dist(SteppableBasePy):
    """
    This class calculates the spatial distribution of polarity with respect to spatial distance.

    """
    def __init__(self,_frequency,_cd):
        SteppableBasePy.__init__(self,_frequency)
        self.cd=_cd
        
    def start(self):   
        # Initializes lists for storing polarity-related metrics.     
        self.PhiB = [] ; self.Phi2B = []
        self.R_matrix = [] ; self.Phi_matrix = []; self.Phibest_matrix = []

        # Initialize PhiB and Phi2B with zero for each bin of size cell diameter (self.cd)
        for R in range(self.cd,self.dim.x,self.cd):
            self.PhiB.append(0); self.Phi2B.append(0)
            
        
    def step(self,mcs):
        # The method iterates through distance bins and computes polarity correlations at each time step
        if (mcs>0):
            L=[]; Lcells=[] #Optimized Algoritm (20x faster)

            # Store center of mass, ID and polarity properties of CYTO cells
            for cell in self.cell_list_by_type(self.CYTO):
                Lcells.append([cell.id,cell.xCOM,cell.yCOM,cell.dict["Ang-x"],cell.dict["Ang-y"],cell.dict["Pol"]])
            
            # Compute pairwise distances between center of masses and store in L
            for id,x0,y0,xx,yy,p in Lcells:
                L.append([]) #L contains a list for every cell
                for id2,x2,y2,xx2,yy2,p2 in Lcells: 
                    if (id!=id2):
                        dx=abs(x0-x2); dy=abs(y0-y2)

                        # Apply periodic boundary conditions
                        if (dx>self.dim.x/2): dx=self.dim.x-dx
                        if (dy>self.dim.y/2): dy=self.dim.y-dy

                        d=sqrt(dx*dx+dy*dy)
                        L[len(L)-1].append([d,xx2,yy2,p2])

                # Sort by distance for efficient nearest-neighbor computations        
                L[len(L)-1].sort()

            # Initialize lists for storing polarity correlations at different distances    
            PhiB=[]; Phi2B=[]; nn=len(L); nC=[0]*nn
            PhiCx=[0]*nn; PhiCy=[0]*nn
            PhiCx2=[0]*nn; PhiCy2=[0]*nn

            # Iterating over the distance bins in the x-direction with steps of size cell diameter.
            for R in range(self.cd,self.dim.x,self.cd):
                Phi=0.0; Phi2=0.0   # Initialize the polarity values for this distance bin
                
                # Loop over each cell to calculate the average polarity at distance R
                for cell in range(nn):
                    n=0     # This will count the number of neighbors within this distance bin
                    
                    # Loop over the neighbors of the current cell (L[cell] contains the neighbors of cell)
                    for d,dx,dy,p in L[cell]:
                        # If the distance between the current cell and its neighbor exceeds R, stop considering further neighbors
                        if (d>R): break 

                        # Increment the count for this neighbor and add the x and y components of the orientation for this cell
                        n+=1; PhiCx[cell]+=dx; PhiCy[cell]+=dy

                        # Accumulate the weighted contributions (based on polarity) for this cell
                        PhiCx2[cell]+=dx*p; PhiCy2[cell]+=dy*p

                     # Keep track of how many neighbors were found within this distance bin    
                    nC[cell]+=n

                    # If there were neighbors within this distance bin, calculate the polarity
                    if (nC[cell]>0):
                        for i in range(n):
                            # Remove processed neighbors for this cell to avoid recalculating
                            del L[cell][0]

                        # Calculate the average polarity (magnitude) for the cell using the accumulated orientation components    
                        Phi+=sqrt(PhiCx[cell]*PhiCx[cell]+PhiCy[cell]*PhiCy[cell])/nC[cell]
                         # Calculate the weighted average polarity 
                        Phi2+=sqrt(PhiCx2[cell]*PhiCx2[cell]+PhiCy2[cell]*PhiCy2[cell])/nC[cell]
                
                self.R_matrix.append(R/self.cd)     # Store the normalized distance
                self.Phi_matrix.append(Phi/nn)      # Store the mean polarity for this distance bin
                
                # Store the mean and weighted polarity for this bin
                PhiB.append(Phi/nn); Phi2B.append(Phi2/nn)

            # After completing the distance bins, update the best observed polarity (PhiB)
            # If the sum of the current PhiB is greater than the previous best PhiB, update PhiB
            if (sum(PhiB)>sum(self.PhiB)):
                self.PhiB = PhiB

            # Similarly, update the best polarity for Phi2 if the current sum of Phi2B exceeds the previous best
            if (sum(PhiB)>sum(self.Phi2B)):
                self.Phi2B = Phi2B

            # Append the best polarity for each distance bin to Phibest_matrix for future reference
            for r in range(len(PhiB)):
                self.Phibest_matrix.append(self.PhiB[r])      
        
    
    def finish(self):
        #Save the computed polarity correlation data to a file ("localpol.txt")
        output_file, path__ = self.open_file_in_simulation_output_folder("localpol.txt", mode='w')
        #output_file.write('R Phi Phi_best\n')
        for i in range(len(self.R_matrix)):
            output_file.write('{} {} {}\n'.format(self.R_matrix[i],self.Phi_matrix[i],self.Phibest_matrix[i]))
               
    
class GraphPol_x_Dist_left_boundary(SteppableBasePy):
    """
    This class calculates the spatial distribution of polarity in relation to distance from the left boundary.

    """
    def __init__(self,_frequency,_cd, _t_relax):
        SteppableBasePy.__init__(self,_frequency)
        self.cd=_cd
        self.t_relax=_t_relax

    def start(self):
        # Initializing arrays to store polarization and angle data
        self.PhiB = [] ; self.Phi2B = []
        self.angle = []
        self.R_matrix = [] ; self.Phi_matrix = []; self.Phibest_matrix = []

        # Setting up bins based on the step size cell diameter for distances along the x-axis
        for R in range(self.cd,self.dim.x,self.cd):
            self.PhiB.append(0); self.Phi2B.append(0)


    def step(self, mcs):
        if mcs > self.t_relax:      # Start analysis after the relaxation time
            # Create bins along the x-axis based on cell diameter step
            bin_edges = list(range(2+self.cd, self.dim.x, self.cd))  
            num_bins = len(bin_edges)

            # Initialize arrays for storing polarization components and cell counts in each bin
            Phi_x_bins = [0] * num_bins
            Phi_y_bins = [0] * num_bins
            cell_counts = [0] * num_bins

            # Assign cells to corresponding bins based on their x-coordinate (COM)
            for cell in self.cell_list_by_type(self.CYTO):
                for bin_idx, bin_edge in enumerate(bin_edges):
                    if cell.xCOM < bin_edge:
                        Phi_x_bins[bin_idx] += cell.dict["Ang-x"]
                        Phi_y_bins[bin_idx] += cell.dict["Ang-y"]
                        cell_counts[bin_idx] += 1


            # Normalize polarization components
            for bin_idx in range(num_bins):
                if cell_counts[bin_idx] > 0:
                    Phi_x_bins[bin_idx] /= cell_counts[bin_idx]
                    Phi_y_bins[bin_idx] /= cell_counts[bin_idx]

            # Initialize an array to store the local polarization magnitude at each bin
            Local_Magnitude = []
            bin_centers = []

            # Calculate the polarization magnitude at each bin
            for bin_idx in range(num_bins):
                magnitude = (Phi_x_bins[bin_idx]**2 + Phi_y_bins[bin_idx]**2)**0.5
                Local_Magnitude.append(magnitude)   # Store the magnitude of polarization for the bin
                bin_centers.append(bin_idx)         # Store the bin index (center)

            # Store the results
            for i in range(len(bin_centers)):
                self.R_matrix.append(bin_centers[i])
                self.Phi_matrix.append(Local_Magnitude[i])

            # Write the results to a file ("localpol.txt") for the current MCS (Monte Carlo Step)
            output_file, path__ = self.open_file_in_simulation_output_folder("localpol"+str(mcs)+".txt", mode='w')

            for i in range(len(self.R_matrix)):
                output_file.write('{} {}\n'.format(self.R_matrix[i],self.Phi_matrix[i]))


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
       
class UniformGrowthSteppable(SteppableBasePy):         
    """
    This class implements cell growth for uniform cell proliferation

    Arguments:
            frequency -- Frequency at which this Steppable runs
            _cd -- cell diameter 
            _tV -- target volume of the cell
            _LamV -- lambda volume of the cell
            _pP -- volume fraction for PROX2 cells
            _pD -- volume fraction for DIST2 cells
            _pC -- volume fraction for CYTO cells
            _growth_rate -- rate of growth for the cells
            t_relax -- relaxation time before growth starts    
    """
    def __init__(self, _frequency, _cd,_tV,_LamV,_pP,_pD,_pC,_growth_rate,t_relax):
        SteppableBasePy.__init__(self, _frequency)
        self.cd=_cd; self.tV=_tV; self.LamV=_LamV; 
        self.pP=_pP; self.pD=_pD; self.pC=_pC; self.pL=1-_pP-_pD-_pC
        self.growth_rate = _growth_rate
        self.t_relax = t_relax
        
    def start(self):
        # Initializes the dictionary to store growth start times for each cluster.
        self.growth_start_times = {}

    def step(self, mcs):
        if(mcs>=self.t_relax):      # Growth starts only after relaxation time has passed                              
            growth_cluster_id_list = []               #list to store cluster ids of cells to grow
               
            for cell in self.cell_list_by_type(self.CYTO):              #goes over all clusters
                if (self.cd<=cell.xCOM<=(self.dim.x - 2*self.cd)):      #check if the cell is within the dimensions of the lattice              
                    if "division_time" in cell.dict:                    #check if the cell was divided before(division time is the sum of the mcs at which the cell was divided 
                        if ((mcs - cell.dict["division_time"]) >= 0):   #and random number between 0 and maximum relaxation time before starting another growth
                            #print("yes")                            
                            growth_cluster_id_list.append(cell.clusterId)       #save the cluster id of the cell to grow only if the mcs is greater than/equal to division time
                            
                    else:
                        growth_cluster_id_list.append(cell.clusterId)   #save the cluster id of the cell to grow if it was not divided before        
                        
                                       
            cells_to_grow = np.unique(growth_cluster_id_list)           #removes duplicate cluster ids 
            
            if (mcs == self.t_relax):                                   #choosing the first cells to grow at relaxation time
                for cluster_to_grow in cells_to_grow:                   #going over all cells to grow
                    time_to_start_grow = randint(self.t_relax, 2*self.t_relax)          #assigning time for each cell to start growing by drawing a random number 
                    self.growth_start_times[cluster_to_grow] = {"time_to_grow":time_to_start_grow}              #storing the growth start times in a dictionary
            
            for cluster_to_grow in cells_to_grow:                       #going over all cells to grow
                if cluster_to_grow in self.growth_start_times:
                    if mcs >= self.growth_start_times[cluster_to_grow]["time_to_grow"]:    #if the current time is greater than the relaxation time                    
                        #estimate the current target volume of the cells
                        cluster_volume = 0          
                        for cell in self.get_cluster_cells(int(cluster_to_grow)):
                            cluster_volume += cell.targetVolume

                        if(cluster_volume<2*self.tV): #if the current target volume is less than twice the volume, increase the volume of the cell
                            #adding volume 
                            cluster_volume+=self.growth_rate
                            S= 4*np.sqrt(cluster_volume)        #estimate the current target surface area of the cell based on its current target volume
                            
                            #assign the target volume to each cell compartment
                            for cell in self.get_cluster_cells(int(cluster_to_grow)):                        
                                if(cell.type==self.CYTO):
                                    cell.targetVolume = (cluster_volume - 1.5*S)
                                elif (cell.type==self.PROX2):
                                    cell.targetVolume = 0.45*S
                                elif (cell.type==self.DIST2):
                                    cell.targetVolume = 0.45*S
                                elif (cell.type==self.LAT):
                                    cell.targetVolume = 0.6*S  
                                    
                else:
                    # If no specific start time is assigned, grow normally
                    cluster_volume = 0          
                    for cell in self.get_cluster_cells(int(cluster_to_grow)):
                        cluster_volume += cell.targetVolume

                    if(cluster_volume<2*self.tV): 
                        cluster_volume+=self.growth_rate    # adding volume
                        S= 4*np.sqrt(cluster_volume)
                        
                        #assign the target volume to each cell compartment
                        for cell in self.get_cluster_cells(int(cluster_to_grow)):                        
                            if(cell.type==self.CYTO):
                                cell.targetVolume = (cluster_volume - 1.4*S)
                            elif (cell.type==self.PROX2):
                                cell.targetVolume = 0.33*1.5*S
                            elif (cell.type==self.DIST2):
                                cell.targetVolume = 0.33*1.5*S
                            elif (cell.type==self.LAT):
                                cell.targetVolume = 0.33*1.5*S  


class MitosisUniformGrowthSteppable(MitosisSteppableClustersBase):
    """ This class implements mitosis for uniform cell proliferation.
    Arguments:
          frequency -- Frequency at which this Steppable runs
            _cd -- cell diameter used for normalizing polarity  
            _tV -- target volume of the cell
            _LamV -- lambda volume of the cell  
            _pP -- proportion of proximal domain
            _pD -- proportion of distal domain
            _pC -- proportion of cytoplasm domain
            _t_relax -- relaxation time for division    
    """
    def __init__(self, frequency,_cd,_tV,_LamV,_pP,_pD,_pC,t_relax):
        MitosisSteppableClustersBase.__init__(self, frequency)
        self.cd=_cd; self.tV=_tV; self.LamV=_LamV; 
        self.pP=_pP; self.pD=_pD; self.pC=_pC; self.pL=1-_pP-_pD-_pC
        self.t_relax = t_relax
        
    def step(self, mcs):                        
        cells_to_divide = []        #  list to store cluster ids of cells to divide
        if(mcs>=self.t_relax):      # Growth starts only after relaxation time has passed
            mitosis_cluster_id_list = []    # List to store cluster ids of cells to divide
            for cell in self.cell_list_by_type(self.CYTO):  
                cluster_cell_list = self.get_cluster_cells(cell.clusterId)
                cluster_volume = 0
                for cell_cmpt in cluster_cell_list:
                    cluster_volume += cell_cmpt.volume

                if(cluster_volume>=2*self.tV): # if the current target volume is greater than twice the volume, divide the cell
                    mitosis_cluster_id_list.append(cell.clusterId) 
                
            for cluster_to_divide in mitosis_cluster_id_list:
                #print(cluster_to_divide,"dividing - yes")    
                self.divide_cluster_along_minor_axis(cluster_to_divide)     # to divide the cell along minor axis
                #self.divide_cluster_random_orientation(cluster_to_divide)   # to divide the cell along random axis
                self.check()                                                # to check if the cell is divided and inherited compartments correctly
                self.check_for_ratios()                                     # to check if the cell is divided and inherited compartments correctly

                for cell_cmpt in self.get_cluster_cells(self.parent_cell.clusterId):    
                    if (cell_cmpt.type==self.CYTO):
                        cell_cmpt.dict["division_time"] = mcs+randint(0,self.t_relax)   # store the division time for the parent cell
                        # print("printing cell id and diviison time for parent cell")
                        # print(cell_cmpt.id, cell_cmpt.dict["division_time"])
                for cell_cmpt in self.get_cluster_cells(self.child_cell.clusterId):
                    if (cell_cmpt.type==self.CYTO):
                        cell_cmpt.dict["division_time"] = mcs+randint(0,self.t_relax)   # store the division time for the child cell
                        # print("printing cell id and diviison time for child cell")
                        # print(cell_cmpt.id, cell_cmpt.dict["division_time"])
 
   
    def check(self):
        for cell in self.cell_list_by_type(self.CYTO):
            cluster_cell_list = self.get_cluster_cells(cell.clusterId)
            cluster_volume = 0
            for cell_cmpt in cluster_cell_list: 
                cluster_volume += cell_cmpt.targetVolume
            if cluster_volume > (2*self.tV + 0.5):                      # check if the target volume of the cell is greater than twice the target volume even after division
                #print("#########################CELL IS WRONG##########################")
                #print("cluster volume",cluster_volume)
                for cell_cmpt in cluster_cell_list:
                    #print(cell_cmpt.type, " volume " ,cell_cmpt.volume," target volume ", cell_cmpt.targetVolume," lambda volume ", cell.lambdaVolume)
                    if(cell_cmpt.type==self.PROX2):
                        cell_cmpt.targetVolume=(self.tV)*self.pP           # assign the target volume to each cell compartment of the newly divided cell which did not update its target volumes 
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
                    #print("%%%%%%new volume%%%%%%")
                    #print(cell_cmpt.type, " volume " ,cell_cmpt.volume," target volume ", cell_cmpt.targetVolume," lambda volume ", cell.lambdaVolume)
        #for cell in self.cells_to_delete:
            #self.delete_cell(cell)
            
    def check_for_ratios(self):                     # to fix the bug that doesn't assign the correct target volume to the newly divided cell
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
            
    
    def update_attributes(self):        # to update the attributes of the newly divided cells 
        self.compartment_list_parent = self.get_cluster_cells(self.parent_cell.clusterId)   # get the list of compartments of the parent cell
        self.compartment_list_child = self.get_cluster_cells(self.child_cell.clusterId)     # get the list of compartments of the child cell
        
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
                    if (new_cell.id==cell_.id):         # check if the cell is one of the compartments of the newly divided cell
                        pList.append([x,y,0])           # save the pixel coordinates

        shuffle(pList)                      # shuffle the list of pixels of the newly divided cell
        cluster_volume = len(pList)         # get the volume of the newly divided cell
        p_frac = self.pP                    # get the volume fraction of the proximal domain
        d_frac = self.pD                    # get the volume fraction of the distal domain
        l_frac = self.pL                    # get the volume fraction of the lateral domain
  
        for x,y,z in pList:                 # randomly assign the cell type to each pixel of the newly divided cell
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
        cytoCell.targetVolume=(self.tV)*(1-p_frac-d_frac-l_frac)            # assign the target volume to each cell compartment
        cytoCell.lambdaVolume=self.LamV*2        
        reassignIdFlag=self.reassign_cluster_id(cytoCell,cell.clusterId)    # reassign the cluster id to the newly divided cell
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


class GrowthOnFrontSteppable(SteppableBasePy):
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
    def __init__(self, _frequency, _cd,_tV,_LamV,_pP,_pD,_pC,_growth_rate,t_relax):
        SteppableBasePy.__init__(self, _frequency)
        self.cd=_cd; self.tV=_tV; self.LamV=_LamV; 
        self.pP=_pP; self.pD=_pD; self.pC=_pC; self.pL=1-_pP-_pD-_pC
        self.growth_rate = _growth_rate
        self.t_relax = t_relax

    def step(self, mcs):
        if(mcs>=self.t_relax):      # Growth starts only after relaxation time has passed
            growth_cluster_id_list = []                     # list to store cluster ids of cells to grow
            shrink_cluster_id_list = []                     # list to store cluster ids of cells to shrink
            for cell in self.cell_list_by_type(self.CYTO):      # goes over all clusters
                if (self.cd<=cell.xCOM<=(self.dim.x - 2*self.cd)): # check if the cell is within the dimensions of the lattice       
                    # calculate the volume of the cluster
                    cluster_volume = 0       
                    for cell_cmpt in self.get_cluster_cells(cell.clusterId):   
                        cluster_volume += cell_cmpt.targetVolume  
                        # check if the cell is on the front (by checking if the cell/compartment has no neighbors)
                        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell_cmpt):
                            if (not neighbor):                
                                growth_cluster_id_list.append(cell.clusterId)
                                break
                    # check if the cell is not on the front and has a volume greater than target volume
                    if (cell.clusterId not in growth_cluster_id_list) and (cluster_volume>self.tV):
                        shrink_cluster_id_list.append(cell.clusterId)   # save the cluster id of the cell to shrink
                                       
            cells_to_grow = np.unique(growth_cluster_id_list)   # removes duplicate cluster ids
            cells_to_shrink = shrink_cluster_id_list            
            
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
                            


class MitosisOnFrontSteppable(MitosisSteppableClustersBase):
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
    def __init__(self, frequency,_cd,_tV,_LamV,_pP,_pD,_pC,t_relax):
        MitosisSteppableClustersBase.__init__(self, frequency)
        self.cd=_cd; self.tV=_tV; self.LamV=_LamV; 
        self.pP=_pP; self.pD=_pD; self.pC=_pC; self.pL=1-_pP-_pD-_pC
        self.t_relax = t_relax
        
    def step(self, mcs):
        cells_to_divide = []
        if(mcs>=self.t_relax):      # Growth starts only after relaxation time has passed
            mitosis_cluster_id_list = []
            for cell in self.cell_list_by_type(self.CYTO):
                cluster_cell_list = self.get_cluster_cells(cell.clusterId)
                cluster_volume = 0
                for cell_cmpt in cluster_cell_list:
                    cluster_volume += cell_cmpt.volume
                if(cluster_volume>=2*self.tV):# :       # if the current target volume is greater than twice the target volume, divide the cell
                    mitosis_cluster_id_list.append(cell.clusterId) 
            for cluster_to_divide in mitosis_cluster_id_list:           # going over all cells to divide
                print(cluster_to_divide,"dividing - yes")    
                #self.divide_cluster_random_orientation(cluster_to_divide)   # to divide the cell along random axis
                self.divide_cluster_along_minor_axis(cluster_to_divide)         # to divide the cell along minor axis
                self.check()        # to check if the cell is divided and inherited compartments correctly
                self.check_for_ratios() # to check if the cell is divided and inherited compartments correctly
    
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
                                    
