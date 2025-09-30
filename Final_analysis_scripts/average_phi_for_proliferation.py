# -*- coding: utf-8 -*-
"""
Created on Tue May  7 09:42:16 2024

@author: abhis
"""

from math import *
import numpy as np
import os
#import matplotlib.pyplot as plt
#plt.rcParams.update({'font.size': 15})

#folder_name_ = "/share/belmonte/athayam/CC3DWorkSpace/Wing_new_final__2025_Feb_10__10_42_23"
#folder_name_ = "/share/belmonte/athayam/CC3DWorkSpace/Wing_new_final__relax_time_10000000__2025_Mar_31__11_30_55"
#folder_name = os.getcwd()
folder_name_ = "/share/belmonte/athayam/CC3DWorkSpace/Wing_new_final__relax_time_10000000__2025_Aug_28__08_53_55"
N = 100

target_values = [40, 100, 200, 300, 400, 500, 600, 700,800,900]
selected_values = []
selected_phi_values = []
average_phi = np.zeros(len(target_values))
n_counter_list = np.zeros(len(target_values))
angle_list = []
for n in range(1,N+1):
    #folder_name = "/share/belmonte/athayam/CC3DWorkSpace/Wing_new_final__2025_Feb_10__10_42_23"
    str_n = "%04d" % n
    folder_name = folder_name_ + "/"+str_n + "/" 
    filename = "globalpol.txt"
    print(folder_name+filename)
    
    if os.path.isfile(folder_name+filename):
        print("yes")
        data = np.loadtxt(folder_name+"globalpol.txt", delimiter=" ")
        time = data[:,0]
        phi = data[:,1]
        n_cells = data[:,2]
        angle = data[:,3]
        n_cells_init = n_cells[0]
        N_cell_list = []
        phi_list = []
        
        for i in range(1,len(n_cells)):
            n_cells_final = n_cells[i]
            if (n_cells_final > n_cells_init):
                N_cell_list.append(n_cells_init)
                phi_list.append(phi[i-1])
            n_cells_init = n_cells_final

        #print(phi_list)
        for j, target_value in enumerate(target_values):
            if N_cell_list:
                closest_value = min(N_cell_list, key=lambda x: abs(target_value-x))
                print(closest_value)
                index = N_cell_list.index(closest_value)
                selected_values.append(N_cell_list[index])
                selected_phi_values.append(phi_list[index])
                average_phi[j] += phi_list[index]
                n_counter_list[j] += 1
        angle_list.append(angle[index])

for i in range(len(target_values)):
    average_phi[i] = average_phi[i]/n_counter_list[i]    
print(target_values,average_phi)

np.savetxt("/share/belmonte/athayam/PCP_REVISION/angles_uniform_t_10000_j_int_all_neg_2.csv",
               angle_list, delimiter=",", header="Final Angle", comments='')




