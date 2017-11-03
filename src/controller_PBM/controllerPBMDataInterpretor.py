#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  23 14:06:49 2017

@author: Chaitanya Sampat
"""
"""
This file reads the data from an initial time step and stores it consequently.
It also compares the last collected data with the intial data and checks if there is a difference
in the properties and instructs whether the PBM should continue, changed to the DEM or that the
system has reached steady state

==================================================================================================
OUTPUT FLAG VALUE STATUS
==================================================================================================
0 - Continue with PBM
1 - Kill PBM and execute DEM with the new files
2 - The system has reached state with the PBM execution, kill all executions
==================================================================================================
"""
import numpy as np
import controllerPBMDataReader as PBM_reader

class controllerPBMDataInterpretor(object):

    def __init__(self, init_ts, compartment, bins1, bins2, PBM_out_path):
        self.initial_ts = np.float(init_ts)
        self.bins1 = bins1
        self.bins2 = bins2
        self.compartments = compartment
        self.PBM_output_path = PBM_out_path
        self.obj_PBM_reader = PBM_reader.controllerPBMDataReader(self.initial_ts, self.compartments, self.bins1, self.bins2, self.PBM_output_path)
        self.initial_d50 = np.array(self.compartments +1)
        self.initial_num_particles = 0.0
        self.init_data_storage(self.initial_ts)
        self.d50_store = np.zeros_like(self.initial_d50)
        self.num_particles = np.zeros((1,2))
        self.new_data_storage(self.initial_ts)
        self.time_index = 0


    # method to calculate the initial d50 and number of particles.
    def init_data_storage(self, init_ts):
        obj_PBM_reader = PBM_reader.controllerPBMDataReader(init_ts, self.compartments, self.bins1, self.bins2, self.PBM_output_path)
        self.initial_d50 = obj_PBM_reader.data_d50_extractor(init_ts)
        self.initial_num_particles = obj_PBM_reader.data_particles_extractor(init_ts)

    # method to store the new d50 and number of particles data from the files being printed
    def new_data_storage(self, init_ts):
        new_ts = np.float(self.obj_PBM_reader.nextfile_time_finder(init_ts))
#        new_d50 = self.initial_d50
#        new_num_particles = 0.0
        temp_new_d50 = np.zeros_like(self.initial_d50)
        temp_new_d50 = self.obj_PBM_reader.data_d50_extractor(new_ts)
        temp_new_num_particles = self.obj_PBM_reader.data_particles_extractor(new_ts)
        self.d50_store = np.append(self.d50_store, temp_new_d50)
        self.num_particles = np.append(self.num_particles, [init_ts, temp_new_num_particles])

    # method to compare the data from the data and decide on the execution
    def data_comparison(self, ts):
        flag = 0
        print(ts)
        temp1 = self.obj_PBM_reader.data_d50_extractor(ts)
        temp2 = self.num_particles[-1]
        particles_check = abs(self.initial_num_particles - temp2) / self.initial_num_particles
        ss_check_time = 10
        # print("Now checking %f with initial time %f"%(ts, self.initial_ts))
        # This checks if the current timestep is greater than the older timestep and whether the properties have changed by 25% so that we kill PBM and start DEM.
        if (ts > self.initial_ts and particles_check < 0.25):
            for i in range(0,self.compartments):
                d50_check = (abs(self.initial_d50[0][i] - temp1[0][i]) / self.initial_d50[0][i])
                if (int(d50_check) > 0.25):
                    flag = 1
                    break
                else:
                    continue
        elif ((ts - self.initial_ts) > ss_check_time):
            flag = 2
        elif ((ts - self.initial_ts) < ss_check_time and particles_check > 0.15):
            flag = 1
        else:
            flag = 0
        return flag

