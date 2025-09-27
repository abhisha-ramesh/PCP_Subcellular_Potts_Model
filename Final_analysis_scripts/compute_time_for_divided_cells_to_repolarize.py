#!/usr/bin/env python3

import os
import numpy as np

# Compute time taken for divided cells to re-polarize after division based on its inner product with P-D axis
def find_first_crossing_from_file(filename1, filename2):
    # Read and parse the file
    data = []
    data_phi = np.loadtxt(filename1, usecols=(0, 1))

    # Separate columns
    time_data = data_phi[:, 0]  # first column
    pol_data = data_phi[:, 1]  # second column
    
    with open(filename1, 'r') as f:
        #print("opening file")
        for line in f:
            # Skip empty lines
            if not line.strip():
                continue
            # Split by whitespace and convert to float/int
            parts = line.split()
            time_val = float(parts[0])
            col2 = float(parts[1])
            col3 = int(parts[2])
            col4 = float(parts[3])
            col5 = float(parts[4])
            data.append([time_val, col2, col3, col4, col5])


    with open(filename2, 'r') as f:
        division_time = []
        for line in f:
            # Skip empty lines
            if not line.strip():
                continue
            division_time.append(int(line.split()[0]))

    start_time = np.mean(division_time)      # find the time of cell division

    threshold = 0.95  # threshold for inner product with P-D axis (0.95 means within ~20 degrees of P-D axis)

    # Find first crossing
    for row in data:
        if row[0] > start_time and row[4] > threshold:  # Check if time is after division and inner product exceeds threshold
            time__ = row[0]-start_time
            return time__

# find the time for 3 different growth rates, 2 different division axes (major and minor) and 3 different criteria
# criteria 1: single cell inside the tissue
# criteria 2: a column of cells inside the tissue
# criteria 03: a column of cells in the front boundary of the tissue

growth_rate = [0.1,0.01,0.001]
axis = [0,1]
criteria = [0,1,2]
for gr in growth_rate:
    for ax in axis:
        for cr in criteria:
            time_list = []
            PHI_list = np.zeros(50000)  #500000 MCS was the max time of simulation (we save the data every 10 MCS)
            n_counter = np.zeros(50000)
            for n in range(1, 101):
                str_n = "%03d" % n
                filename1 = f'/share/belmonte/athayam/CC3DWorkSpace/cell_non_auto_time_for_dividing_cells_to_polarize_along_P_D/__axis_{str(ax)}__gr_{str(gr)}__gr_criteria_{str(cr)}/0{str_n}/globalpol.txt'
                filename2 = f'/share/belmonte/athayam/CC3DWorkSpace/cell_non_auto_time_for_dividing_cells_to_polarize_along_P_D/__axis_{str(ax)}__gr_{str(gr)}__gr_criteria_{str(cr)}/0{str_n}/time_of_division.txt'
                if os.path.isfile(filename1):
                    if os.path.isfile(filename2):
                        result = find_first_crossing_from_file(filename1, filename2)
                        time_list.append(result)
            print("axis_"+str(ax)+"_criteria_"+str(cr)+"_growth_"+str(gr)+"=", time_list)
            
