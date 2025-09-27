#!/usr/bin/env python3

import numpy as np
import os
import csv
import sys



def T_ave(N_sim,N):
    n_files = 0
    Z_sum = np.zeros(99000)
    n_sum = 0
    Time = np.zeros(99000)
    N_cells_sum = np.zeros(99000)
    time_for_50 = []
    for n in range(1,N_sim+1):
        str_n = "%03d" % n
        #filename = '/share/belmonte/athayam/CC3DWorkSpace/PCP__2025_Aug_04__19_24_20/__Nx_'+str(N) + '/0' + str_n + '/globalpol.txt'i
        filename = f'/share/belmonte/athayam/CC3DWorkSpace/system_size_90x10_lb_J_pd_neg_1/__N_{str(N)}/0{str_n}/globalpol.txt'
        #filename = f'/share/belmonte/athayam/CC3DWorkSpace/periodic_system_size_data_with_adhesion_fix/PCP__N_{str(N)}__2025_Feb_15__22_14_16/0{str_n}/globalpol.txt'
        if os.path.isfile(filename):
            data = np.loadtxt(filename, delimiter=' ')
            #print(filename)
            X = data[:, 0]             
            Z = data[:, 1]
            N_cells = data[:, 2]
            for i in range(len(X)):
                if Z[i] > 0.5:
                    time_for_50.append(X[i])
                    break
        else:
            n_files += 1
    print("Failed simulations",n_files)
    for i in range(len(X)):
        if len(time_for_50)>0:
            ave_time_for_50 = np.mean(time_for_50)
            std_time = np.std(time_for_50)
            se_time = np.std(time_for_50)/np.sqrt(len(time_for_50))
        else:
            ave_time_for_50 = 0
            std_time = 0
            se_time = 0
        
    return ave_time_for_50, std_time, se_time, len(time_for_50)

if __name__ == "__main__":
    N_list = [10,20,30,50,70,90]
    av_t = []
    std_t = []
    se_t = []
    num_sims = []
    for N in N_list:
        ave_time, std_time, se_time, n_sims = T_ave(100,N)
        av_t.append(ave_time)
        std_t.append(std_time)
        se_t.append(se_time)
        num_sims.append(n_sims)
    np.savetxt("/share/belmonte/athayam/PCP_REVISION/ave_time_lb_Nx_10_j_pd_neg_1.csv", np.c_[N_list, av_t, std_t, se_t, num_sims], delimiter=",",
               header="N, Average Time for 50, STD Time, SE Time, Number of Simulations")

