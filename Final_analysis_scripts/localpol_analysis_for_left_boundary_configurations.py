#!/usr/bin/env python3

import numpy as np
import os


# Compute average local polarization as a function of distance from left boundary for different system sizes
def local_pol(N_sim,n_cells):
    local_pol_data = np.zeros((N_sim, n_cells+n_cells))
    cell_number_list = np.zeros(n_cells+n_cells)
    n_files = 0
    str_n_cells = str(n_cells)
    for n in range(1, N_sim + 1):
        str_n = "%03d" % n
        filename = f'/share/belmonte/athayam/CC3DWorkSpace/PCP__2025_Sep_04__20_52_34/__N_{str(n_cells)}/0{str_n}/localpol.txt'
        # filename = f'0{str_n}/localpol.txt'
        
        if os.path.isfile(filename):
            # print("yes")
            data = np.loadtxt(filename, delimiter=' ')
            n_files += 1
            if data.size > 0:
                for n_cell in range(1,n_cells+n_cells):
                    local_pol_data[n-1,(n_cells+n_cells - n_cell-1)] = data[-n_cell,1]
                    cell_number_list[n_cell-1] = n_cell
                    
    local_pol_sum = np.sum(local_pol_data, axis=0, keepdims=True)
    
    if (n_files != 0):    
        localpol_data_final = (local_pol_sum/n_files)        
        L = localpol_data_final.flatten()
        # print(cell_number_list)
        # print(L)
    # Save final angles to a CSV file
    np.savetxt("/share/belmonte/athayam/PCP_REVISION/localpol_left_boundary_N_"+str(n_cells)+".csv", np.c_[cell_number_list, L], delimiter=",", header="radius, localpol")
    # np.savetxt("localpol_periodic_N_30.csv", 
    #            np.c_[cell_number_list, L], delimiter=",", header="radius, localpol")
    
if __name__ == "__main__":
    system_size_list = [4,8,12,16,20,24,30]
    for N in system_size_list:
        local_pol(100, N)  # Using 300 simulations



