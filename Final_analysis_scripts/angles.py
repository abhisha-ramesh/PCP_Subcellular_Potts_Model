#!/usr/bin/env python3

import numpy as np
import os


# Extract final angle of polarization from multiple simulations and save to CSV
def T_final_angles(N_sim):
    final_angles = []  # List to store final angles

    for n in range(1, N_sim + 1):
        str_n = "%03d" % n
        filename = f'/share/belmonte/athayam/CC3DWorkSpace/PCP__2025_Sep_03__22_08_26/0{str_n}/globalpol.txt'

        if os.path.isfile(filename):
            data = np.loadtxt(filename, delimiter=' ')
            # Check if data is not empty and take the last angle
            print(filename)
            if data.size > 0:
                if (len(data[:,1]) == 199000):
                    #if data[-1,1]>=0.9:
                    final_angle = data[-1, 3]  #  angle is in the fourth column
                    final_angles.append(final_angle)
                
    # Save final angles to a CSV file
    output_path = "/share/belmonte/athayam/PCP_REVISION/angles_periodic_8x8.csv"
    np.savetxt(output_path, 
               final_angles, delimiter=",", header="Final Angle", comments='')

if __name__ == "__main__":
    T_final_angles(100)  # Using 100 simulations

