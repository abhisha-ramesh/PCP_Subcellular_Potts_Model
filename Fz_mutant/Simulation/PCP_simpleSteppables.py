from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
from math import *
from random import *
import time
from pathlib import Path

# Initial Condition Steppable with Left Boundary Signal and Periodic Boundary along y-direction
class InitCondY_bXl(SteppableBasePy):
    def __init__(self, _frequency, _cd, _tV, _LamV, _pP, _pD, _pC):
        SteppableBasePy.__init__(self, _frequency)
        self.cd = _cd
        self.tV = _tV
        self.LamV = _LamV
        self.pP = _pP
        self.pD = _pD
        self.pC = _pC
        self.pL = 1 - _pP - _pD - _pC

    def start(self):
        # Initialization of cells with a specific signal on the left boundary
        signal = self.new_cell(self.DIST)
        self.cellField[0, 0, 0] = signal
        for y in range(self.dim.y):
            x = 0
            self.cellField[x, y, 0] = signal
            x = 1
            self.cellField[x, y, 0] = signal
        signal.targetVolume = signal.volume
        signal.lambdaVolume = 100000

        # Drawing cells 
        yy = self.cd / 2
        for x0 in range(2, self.dim.x - self.cd, self.cd):
            yy += self.cd / 2
            for y0 in range(0, self.dim.y, self.cd):
                x = x0
                y = (yy + y0) % self.dim.y
                cell = self.new_cell(self.CYTO)
                cell.targetVolume = self.tV
                cell.lambdaVolume = self.LamV
                for x in range(self.cd):
                    x = x0 + x
                    for y in range(self.cd):
                        y = (yy + y0 + y) % self.dim.y
                        self.cellField[x, y, 0] = cell

        # Setting cell properties
        for cell in self.cell_list_by_type(self.CYTO):
            cell.targetVolume = self.tV
            cell.lambdaVolume = self.LamV

        self.get_all_pixels()
        self.P = 0
        self.D = 0
        for cell in self.cell_list_by_type(self.CYTO):
            self.pixelate(cell)

    # Pixelation method to uniformly distribute subcellular compartments inside cells
    def pixelate(self, cell):
        pCell = self.new_cell(self.PROX)
        dCell = self.new_cell(self.DIST)
        lCell = self.new_cell(self.LAT)
        pList = []

        for pixel in self.L[cell.id]:
            pList.append(pixel)
        shuffle(pList)
        for x, y, z in pList:
            r = random()
            if r < self.pP:
                self.cell_field[x, y, z] = pCell
            elif r < self.pP + self.pD:
                self.cell_field[x, y, z] = dCell
            elif r < self.pP + self.pD + self.pL:
                self.cell_field[x, y, z] = lCell

        cell.targetVolume = self.tV * self.pC
        cell.lambdaVolume = self.LamV * 2
        pCell.targetVolume = self.tV * self.pP
        pCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(pCell, cell.clusterId)
        dCell.targetVolume = self.tV * self.pD
        dCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(dCell, cell.clusterId)
        lCell.targetVolume = self.tV * self.pL
        lCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(lCell, cell.clusterId)

        self.P += pCell.volume
        self.D += dCell.volume

    def get_all_pixels(self):
        self.L = {}
        for cell in self.cell_list:
            self.L[cell.id] = []
        for x, y, z in self.every_pixel():
            cell = self.cell_field[x, y, 0]
            if cell:
                self.L[cell.id].append([x, y, 0])

# Metrics to calculate and store global cell polarization
class CellPol(SteppableBasePy):
    def __init__(self, _frequency, _cd):
        SteppableBasePy.__init__(self, _frequency)
        self.cd = _cd

    def start(self):
        for cell in self.cell_list_by_type(self.CYTO):
            cell.dict["Ang-x"] = 0
            cell.dict["Ang-y"] = 0

    def step(self, mcs):
        for cell in self.cell_list_by_type(self.CYTO):
            compList = self.get_cluster_cells(cell.clusterId)
            xP = 0.0
            yP = 0.0
            xD = 0.0
            yD = 0.0
            ok = 0
            p = 0
            for cell2 in compList:
                if cell2.type == self.PROX:  # proximal
                    xP = cell2.xCOM
                    yP = cell2.yCOM
                    ok += 1
                elif cell2.type == self.DIST:  # distal
                    xD = cell2.xCOM
                    yD = cell2.yCOM
                    ok += 1
                    if cell2.type == self.DIST:
                        p += 1
            if ok == 2:
                dx = xP - xD
                dy = yP - yD
                if abs(dy) > self.dim.y / 2:
                    dy = -(self.dim.y - abs(dy)) * dy / abs(dy)
                if abs(dx) > self.dim.x / 2:
                    dx = -(self.dim.x - abs(dx)) * dx / abs(dx)
                d = sqrt(dx * dx + dy * dy)
                if d > 0:
                    cell.dict["Ang-x"] = dx / d
                    cell.dict["Ang-y"] = dy / d
                    cell.dict["Pol"] = d / self.cd
                else:
                    cell.dict["Ang-x"] = 0
                    cell.dict["Ang-y"] = 0
                    cell.dict["Pol"] = 0
            else:
                cell.dict["Ang-x"] = 0
                cell.dict["Ang-y"] = 0
                cell.dict["Pol"] = 0

# Steppable to plot global polarization over MCS
class GraphGlobalPol(SteppableBasePy):
    def __init__(self, _frequency, _cd):
        SteppableBasePy.__init__(self, _frequency)
        self.cd = _cd

    def start(self):
        # Initialization for global polarization plot
        self.M = []
        self.Phi_matrix = []
        self.n_cell_matrix = []
        for cell in self.cell_list_by_type(self.CYTO):
            cell.dict["Ang-x"] = 0
            cell.dict["Ang-y"] = 0

        self.plot_win = self.add_new_plot_window(
            title='Global Polarization',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Phi',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=False
        )

        self.plot_win.add_plot("Pol", style='Lines', color='red', size=5)

    def step(self, mcs):
        if mcs >= 0:
            n = 0
            x = 0.0
            y = 0.0
            for cell in self.cell_list_by_type(self.CYTO):
                n += 1
                x += cell.dict['Ang-x']
                y += cell.dict['Ang-y']
            Phi = 1
            Phi2 = 1
            if n > 0:
                Phi = sqrt(x * x + y * y) / n  # average of the normalized distance between prox-dist compartments
            self.M.append(mcs)
            self.Phi_matrix.append(Phi)

            # Adding data points to the plot
            self.plot_win.add_data_point("Pol", mcs, Phi)

    def finish(self):
        # Save global polarization data to a file at the end of simulation
        output_file, path__ = self.open_file_in_simulation_output_folder("globalpol.txt", mode='w')
        for i in range(len(self.M)):
            output_file.write('{} {} \n'.format(self.M[i], self.Phi_matrix[i]))

# Steppable to create an orientation field vector based on cell polarization
class OrientationField(SteppableBasePy):
    def __init__(self, _frequency=1):
        SteppableBasePy.__init__(self, _frequency)
        self.create_vector_field_cell_level_py("Orientation")

    def step(self, mcs):
        # Creating orientation field based on cell polarization
        field = self.field.Orientation.clear()
        for cell in self.cell_list_by_type(self.CYTO):
            if cell.dict["Pol"] != 0:
                field = self.field.Orientation
                field[cell] = [-cell.dict["Ang-x"], -cell.dict["Ang-y"], 0]


class VectorField(SteppableBasePy):
    def __init__(self,_frequency=1):
        SteppableBasePy.__init__(self,_frequency)
        
    def step(self,mcs):
        self.X = []; self.Y = []; self.X_com = []; self.Y_com = []; 
        for cell in self.cell_list_by_type(self.CYTO):
            if (cell.dict["Pol"]!=0):
                x = cell.dict['Ang-x']; y = cell.dict['Ang-y']
                self.Y.append(-y)
                self.X.append(-x)
                self.Y_com.append(cell.yCOM)
                self.X_com.append(cell.xCOM)
                
        output_file, path__ = self.open_file_in_simulation_output_folder("vector_field_data_"+str(mcs)+".txt", mode='w')
        for i in range(len(self.X)):
            output_file.write('{} {} {} {}\n'.format(self.X_com[i],self.Y_com[i],self.X[i],self.Y[i]))                       
                
class MutantPhenotype(SteppableBasePy):
    def __init__(self, _frequency, _cd, _tV, _LamV, _pP, _pD, _pC):
        SteppableBasePy.__init__(self, _frequency)
        # Parameters for cell initialization
        self.cd = _cd
        self.tV = _tV
        self.LamV = _LamV
        self.pP = _pP
        self.pD = _pD
        self.pC = _pC
        self.pL = 1 - _pP - _pD - _pC
        
    def step(self, mcs):
        if(mcs==50000):
            cell_list = []
            for cell in self.cell_list_by_type(self.CYTO):
                cell_list.append(cell.id)
            #mut_cell_id = choice(cell_list)
            mut_cell_id_list = [34,35,36,28,27,42,43]
            for mut_cell_id in mut_cell_id_list:
                self.create_mut_cell(mut_cell_id)
            cell = self.fetch_cell_by_id(1)
            self.delete_cell(cell)
            
    def create_mut_cell(self,mut_cell_id):
        cluster_cell_list = self.get_cluster_cells(mut_cell_id)
        pCell = self.new_cell(self.MUT)
        pList = []
        self.get_all_pixels()
        for cell_cmpt in cluster_cell_list:
            if(cell_cmpt.type==self.PROX):
                for pixel in self.L[cell_cmpt.id]:
                    pList.append(pixel)
            elif(cell_cmpt.type==self.DIST):
                for pixel in self.L[cell_cmpt.id]:
                    pList.append(pixel)
            elif(cell_cmpt.type==self.LAT):
                for pixel in self.L[cell_cmpt.id]:
                    pList.append(pixel)
        for x, y, z in pList:
            self.cell_field[x, y, z] = pCell
        pCell.targetVolume = self.tV * (1-self.pC)
        pCell.lambdaVolume = self.LamV * 2
        reassignIdFlag = self.reassign_cluster_id(pCell, mut_cell_id)        
        
    def get_all_pixels(self):
        # Method to get all pixels in the simulation without using PixelTracker Plugin
        self.L = {}
        for cell in self.cell_list:
            self.L[cell.id] = []
        for x, y, z in self.every_pixel():
            cell = self.cell_field[x, y, 0]
            if cell:
                self.L[cell.id].append([x, y, 0])    
        

class piff_generator(SteppableBasePy):
    def __init__(self,_frequency, time, _t_piff):
        SteppableBasePy.__init__(self,_frequency)
        self.time=time; self.t_piff=_t_piff
        
    def step(self,mcs):
        self.get_all_pixels()               #get pixels without using pixel tracker
        self.save_piff(self.t_piff+mcs)                #printing piff files 
        
    def finish(self):
        self.get_all_pixels()               #get pixels without using pixel tracker
        self.save_piff(self.time+self.t_piff)                #printing piff files 
                    
                    
    def save_piff(self, mcs):
        out_folder = self.output_dir
        FileName = out_folder+"/PiffFile_"+str(mcs)+".piff"
        piffPath=Path(out_folder).joinpath(FileName)
        with open(piffPath, 'a') as fout:
            fout.write("Include Clusters \n")
            for i in self.L:
                pixel_list = self.L[i]                
                cell = self.fetch_cell_by_id(i)   
                name = self.get_type_name_by_cell(cell)
                for pixel in pixel_list:
                    x = pixel[0]; y = pixel[1]; z = pixel[2]          
                    fout.write("%d %d %s %d %d %d %d %d %d \n" % (cell.clusterId, cell.id, name, x, x, y, y, z, z))
    
    def get_all_pixels(self):
        self.L = {}   
        for cell in self.cell_list:
            self.L[cell.id] = []
        for x, y, z in self.every_pixel():
            cell = self.cell_field[x, y, 0]       
            if cell:
                self.L[cell.id].append([x,y,0])                
        
        
        
        
        
        
        