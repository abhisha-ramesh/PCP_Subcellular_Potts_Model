
#!/usr/bin/env python3

import numpy as np
import os
import csv
import sys

# Compute average and standard error of global polarization (phi) over multiple simulations (for systems without proliferation)

def T_ave(N_sim):
    n_files = 0
    datapoints = 999000
    Phi_sum = np.zeros(datapoints)
    Phi_sq_sum = np.zeros(datapoints)
    n_sum = np.zeros(datapoints)
    Time = np.zeros(datapoints)
    Phi_std = np.zeros(datapoints)
    Phi_SE = np.zeros(datapoints)
    for n in range(1,N_sim+1):
        str_n = "%03d" % n
        filename = '/share/belmonte/athayam/CC3DWorkSpace/PCP__2025_Sep_05__17_10_08' + '/0' + str_n + '/globalpol.txt'
        if os.path.isfile(filename):
            data = np.loadtxt(filename, delimiter=' ')
            #print(filename)
            Time_data = data[:, 0]             
            Phi_data = data[:, 1]

            if (len(Time_data) >= datapoints):
                for i in range(datapoints):
                    Phi_sum[i] += Phi_data[i]
                    Phi_sq_sum[i] += Phi_data[i]**2
                    n_sum[i] += 1
                    Time[i] = Time_data[i]
        else:
            n_files += 1
    print("Failed simulations",n_files)
    for i in range(datapoints):
        if n_sum[i]>0:
            Phi_sum[i] = Phi_sum[i]/n_sum[i]
            Phi_std[i] = np.sqrt((Phi_sq_sum[i]/n_sum[i]) - (Phi_sum[i]**2))
            Phi_SE[i] = Phi_std[i]/np.sqrt(n_sum[i])
        
    return Time, Phi_sum, Phi_SE

if __name__ == "__main__":
    X, Z_sum, Z_std = T_ave(100)
    output_path = "/share/belmonte/athayam/PCP_REVISION/lb_8x8_new_increased_time.csv"
    np.savetxt(output_path, np.c_[X, Z_sum, Z_std], delimiter=",",
               header="Time, Average Phi, Std, SE")

