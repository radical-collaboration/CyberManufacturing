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
    def __init__(self, prev_timestamp, compartments, bins1, bins2):
        self.p_ts = prev_timestamp # saves the last timestamp from the PBM
        self.compartments= compartments # number of comparments inside the PBM
        self.bins1 = bins1 # number of segregations of bins of solid 1 in PBM
        self.bins2 = bins2 # number of segregations of bins of solid 2 in PBM
        self.d50_overtime = np.array((1,compartments))
        self.last_timeindex = 2


    # function to look for the next file in the PBM output folder after the PBM time step
    def nextfile_time_finder(self, last_ts):
        output_PBM_path = os.getcwd()
        output_PBM_path = os.chdir('../twoway_PBM/csvDump/')
        output_PBM_path = os.getcwd()
        file_list_particles = glob.glob("particles_*.csv")
        file_list_particles.sort()
        file_list_d50 = glob.glob("d50_*.csv")
        file_list_d50.sort()
        new_ts = 0.0
        for file1 in file_list_d50:
            try:
                d50_nextfile_timestamp = re.search('d50_(.+?).csv', file1).group(1)
                particles_nextfile = 'particles_' + d50_nextfile_timestamp + '.csv'
                d50_nextfile = 'd50_' + d50_nextfile_timestamp + '.csv'
                if (os.path.isfile(particles_nextfile) and os.path.isfile(d50_nextfile) and np.float(d50_nextfile_timestamp) > last_ts):
                    new_ts = np.float(d50_nextfile_timestamp)
                    break
                else:
                    continue
            except AttributeError:
                d50_nextfile = ''
                particles_nextfile = ''
                #time.sleep(0.2)
                '''
            try:
                particles_nextfile = re.search('particles_(+.?).csv', file1).group(1)
            except:
                #time.sleep(0.2)
                particles_nextfile = ''
                '''
            if (file1 != file_list_d50[-1]):
                continue
            else:
                time.sleep(0.1)
        return new_ts

    # method to extract data from the d50 csv file and save the data uptill a given time step
    def data_d50_extractor(self, last_ts):
        new_ts = self.nextfile_time_finder(last_ts)
        filetoread_d50 = 'd50_' + str(new_ts) + '.csv'
        temp_d50 = np.zeros(1,self.compartments)
        with open(filetoread_d50, 'rb') as d50_current_file:
            read_d50 = pd.read_csv(d50_current_file, header = 0)
            ts_len = len(read_d50.Time)
            if self.timeindex + 1
            for x in range(1, ts_len):
                for i in range(0,self.compartments):
                    if np.isnan(read_d50.iloc[x][i]):
                        temp_d50[i] = 0
                    else:
                        temp_d50[i] = read_d50.iloc[x][i]
                    self.d50_overtime.append(temp_d50)











a = controllerPBMDataReader(10,4,16,16)
t = a.nextfile_time_finder(10.08)

