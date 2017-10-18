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

# class to locate the PBM dump file, identify the timestamp and store the values of the same

class controllerPBMDataReader(object):

    #Constructor
    def __init__(self, prev_timestamp, compartments, bins1, bins2):
        self.p_ts = prev_timestep # saves the last timestamp from the PBM
        self.compartments= compartments # number of comparments inside the PBM
        self.bins1 = bins1 # number of segregations of bins of solid 1 in PBM
        self.bins2 = bins2 # number of segregations of bins of solid 2 in PBM


    # function to look for the next file in the PBM output folder after the PBM time step
    def file_finder(self, last_ts):
        output_PBM_path = os.getcwd()
        output_PBM_path = os.chdir('../twoway_PBM/csvDump/')
        output_PBM_path = os.getcwd()
        file_list = glob.glob("*.csv")
        for file1 in file_list:
            try:
                d50_find_nextfile = re.search('d50_(.+?).csv', file1).group(1)
            except AttributeError:
                d50_nextfile = ''
                #time.sleep(0.2)
            try:
                particles_nextfile = re.search('particles_(+.?).csv', file1).group(1)
            except:
                #time.sleep(0.2)
                particles_nextfile = ''
            if (np.float(d50_nextfile) > last_ts and np.float(particle_nextfile) > last_ts):
                break
            elif (file1 != file_list[-1]):
                continue
            else:
                time.sleep(0.1)

