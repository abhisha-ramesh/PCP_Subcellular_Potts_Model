# -*- coding: utf-8 -*-
"""
Created on Tue May  7 09:42:16 2024

@author: abhis
"""

from math import *
import numpy as np
import os

# Find average phi at specific cell counts across multiple simulations with proliferation
# Also extract corresponding mean tissue angles of polarization

folder_name_ = "/share/belmonte/athayam/CC3DWorkSpace/Wing_new_final__2025_Jul_30__11_13_21/__relax_time_10000000"
N = 100 #number of simulations  

# Target cell counts for analysis
target_values = [40,100,200,300,400,500,600,700,800,900]

# Lists to store selected values and their corresponding phi values
selected_values = []
selected_phi_values = []

# Initialize arrays for averaging
average_phi = np.zeros(len(target_values))
n_counter_list = np.zeros(len(target_values))
angle_data_list = np.zeros((N,len(target_values)))

# Lists to store final cell counts and phi values from each simulation
n_count_list = np.zeros(2000000)
n_count_plus_list = np.zeros(2000000)
n_count_minus_list = np.zeros(2000000)

average_data_list = np.zeros((2000000,5))


for n in range(1,N+1):
    str_n = "%04d" % n
    folder_name = folder_name_ + "/"+str_n + "/" 
    filename = "globalpol.txt"
    # print(folder_name+filename)
    
    if os.path.isfile(folder_name+filename):
        # print("yes")
        data = np.loadtxt(folder_name+"globalpol.txt", delimiter=" ")   #load data from file
        time = data[:,0]
        phi = data[:,1]
        n_cells = data[:,2]
        angles = data[:,3]

        n_cells_init = n_cells[0]       #initial number of cells
        N_cell_list = []
        phi_list = []
        angle_list = []
        
        # Accumulate data for averaging
        for i in range(len(n_cells)):
            average_data_list[i,0] += time[i]
            average_data_list[i,1] += phi[i]
            average_data_list[i,2] += n_cells[i]
            n_count_list[i] += 1
            if angles[i]>=0:                            #separate counts for positive and negative angles
                average_data_list[i,3] += (angles[i])
                n_count_plus_list[i] += 1
            else:
                average_data_list[i,4] += (angles[i])
                n_count_minus_list[i] += 1

        # Find instances where cell count increases and record corresponding phi values
        for i in range(1,len(n_cells)):
            n_cells_final = n_cells[i]
            if (n_cells_final > n_cells_init):
                N_cell_list.append(n_cells_init)
                phi_list.append(phi[i-1])
            n_cells_init = n_cells_final

        # For each target cell count, find the closest actual cell count and record the corresponding phi value
        for j, target_value in enumerate(target_values):
            if N_cell_list:
                closest_value = min(N_cell_list, key=lambda x: abs(target_value-x))
                print(closest_value)
                index = N_cell_list.index(closest_value)
                selected_values.append(N_cell_list[index])
                selected_phi_values.append(phi_list[index])
                average_phi[j] += phi_list[index]
                n_counter_list[j] += 1


# find average phi at each target cell count
for i in range(len(target_values)):
    average_phi[i] = average_phi[i]/n_counter_list[i]    
print(target_values,average_phi)    #average phi at each target cell count





