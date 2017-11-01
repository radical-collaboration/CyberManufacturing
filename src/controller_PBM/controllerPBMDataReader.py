#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  17 16:19:29 2017

@author: Chaitanya Sampat
"""
'''
This file accesses the folder that contains the output files of the PBM and waits for the data
file to be dumped from the PBM execution and reads the d50 and # of particle data and provides the
information to the interpretor class for further processing
'''
import numpy as np
import os
import glob
import re
import time
import pandas as pd

# class to locate the PBM dump file, identify the timestamp and store the values of the same

class controllerPBMDataReader(object):

    #Constructor
    def __init__(self, prev_timestamp, compartments, bins1, bins2, PBM_out_path):
        self.p_ts = prev_timestamp # saves the last timestamp from the PBM
        self.compartments= compartments # number of comparments inside the PBM
        self.bins1 = bins1 # number of segregations of bins of solid 1 in PBM
        self.bins2 = bins2 # number of segregations of bins of solid 2 in PBM
        self.d50_overtime = np.zeros((1,compartments +1))
        self.last_timeindex = 2
        self.no_particles_overtime = np.zeros(1)
        self.PBM_output_path = PBM_out_path



    # function to look for the next file in the PBM output folder after the PBM time step
    def nextfile_time_finder(self, last_ts):
        particles_file_list = glob.glob(self.PBM_output_path + "particles_*.csv")
        num_start_particles = re.search("particles_",particles_file_list[0]).end()
        num_end_particles = re.search("0.csv",particles_file_list[0]).start()
        # the generated file list is sorted by extracting the time step from the file name and then sorted according to the time step
        if (num_end_particles - num_start_particles == 8):
            particles_file_list = sorted(particles_file_list, key=lambda x: float(x[num_start_particles:num_end_particles]))
        else:
            particles_file_list = sorted(particles_file_list, key=lambda x: float(x[num_start_particles:(num_end_particles - 1)]))
        d50_file_list = glob.glob(self.PBM_output_path + "d50_*.csv")
        # file_list_d50 = sorted(file_list_d50, key=lambda x: float(x[-13:-5]))
        num_end_d50 = re.search("0.csv",d50_file_list[0]).start()
        num_start_d50 = re.search("d50_",d50_file_list[0]).end()
        if (num_end_d50 - num_start_d50 == 8):
            d50_file_list = sorted(d50_file_list, key=lambda x: float(x[num_start_d50:num_end_d50]))
        else:
            d50_file_list = sorted(d50_file_list, key=lambda x: float(x[num_start_d50:(num_end_d50 - 1)]))
        # file_list_particles = sorted(file_list_particles, key=lambda x: float(x[-13:-5]))
#        file_list_particles.sort()
        
#        file_list_d50.sort()
        new_ts = 0
        for x,file1 in enumerate(d50_file_list):
            d50_nextfile_timestamp = re.search(self.PBM_output_path + 'd50_(.+?).csv', file1).group(1)
            # print(d50_nextfile_timestamp)
            particles_nextfile = self.PBM_output_path + 'particles_' + d50_nextfile_timestamp + '.csv'
            d50_nextfile = self.PBM_output_path + 'd50_' + d50_nextfile_timestamp + '.csv'
            if (os.path.isfile(particles_nextfile) and os.path.isfile(d50_nextfile) and np.float(d50_nextfile_timestamp) > last_ts):
                new_ts = np.float(d50_nextfile_timestamp)
                break
            elif (file1 == d50_file_list[-1]):
                # once the iteration reaches the end of the list the program waits for half a second and then looks for newer generated files inside the folder and generates a new list
                print("Waiting for the file to be printed")
                time.sleep(0.6)
                d50_file_list = []
                d50_file_list = glob.glob(self.PBM_output_path + "d50_*.csv")
                # file_list_d50 = sorted(file_list_d50, key=lambda x: float(x[-13:-5]))
                num_end_d50 = re.search("0.csv",d50_file_list[0]).start()
                num_start_d50 = re.search("d50_",d50_file_list[0]).end()
                if (num_end_d50 - num_start_d50 == 8):
                    d50_file_list = sorted(d50_file_list, key=lambda x: float(x[num_start_d50:num_end_d50]))
                else:
                    d50_file_list = sorted(d50_file_list, key=lambda x: float(x[num_start_d50:(num_end_d50 - 1)]))
            else:
                continue                
        return new_ts

    # method to extract data from the d50 csv file and save the data uptill a given time step
    def data_d50_extractor(self, curr_ts):
        while(curr_ts == 0):
            curr_ts = self.nextfile_time_finder(curr_ts)
        curr_ts = "%f" %curr_ts
        filetoread_d50 = self.PBM_output_path + 'd50_' + curr_ts + '.csv'
        temp_d50 = np.zeros((1,self.compartments + 1))

        with open(filetoread_d50, 'rb') as d50_current_file:
            read_d50 = pd.read_csv(d50_current_file, header = 0)
            ts_len = len(read_d50.Time)
            for x in range(0, ts_len):
                for i in range(0,self.compartments + 1):
                    if np.isnan(read_d50.iloc[x][i+1]):
                        temp_d50[0][i] = 0
                    else:
                       temp_d50[0][i] = read_d50.iloc[x][i+1]
                self.d50_overtime = np.vstack([self.d50_overtime, temp_d50])
        return temp_d50


    # method to extract data from the particles csv file and save the data uptill a given time step
    def data_particles_extractor(self, curr_ts):
        if(curr_ts == 0):
            curr_ts = self.nextfile_time_finder(curr_ts)
        curr_ts = "%f" %curr_ts
        filetoread_particles = self.PBM_output_path + 'particles_' + curr_ts + '.csv'
        total_num_particles = 0
        with open(filetoread_particles, 'rb') as particles_current_file:
            read_particles = pd.read_csv(particles_current_file, header = 0)
            ts_len = len(read_particles.value)
            for x in range(1, ts_len):
                total_num_particles += read_particles.iloc[x][3]
            self.no_particles_overtime = np.append(self.no_particles_overtime, total_num_particles)
        return total_num_particles
