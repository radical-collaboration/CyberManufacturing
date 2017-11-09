#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  30 13:41:27 2017

@author: Chaitanya Sampat

-------------------------------------------------------------------------------------------------
This method copies the dummy test files from the copypath to the 

"""

import numpy as np
import sys
import random
import time 
import os
import glob
import shutil
import re

def dummy_PBM(pasteTo_path, last_timestep_index, test_number):
    d50_file_list = glob.glob(os.getcwd() + '/dummy_DEM_PBM' + '/test_%d/d50_*.csv'%test_number)
    particles_file_list = glob.glob(os.getcwd() + '/dummy_DEM_PBM' + '/test_%d/particles_*.csv'%test_number)
    # print(d50_file_list)
    # print(particles_file_list)
    num_start_d50 = re.search("d50_",d50_file_list[0]).end()
    num_start_particles = re.search("particles_",particles_file_list[0]).end()
    num_end_d50 = re.search(".csv",d50_file_list[0]).start()
    num_end_particles = re.search(".csv",particles_file_list[0]).start()
    if (num_end_d50 - num_start_d50 == 8):
        d50_file_list = sorted(d50_file_list, key=lambda x: float(x[num_start_d50:num_end_d50]))
    else:
        d50_file_list = sorted(d50_file_list, key=lambda x: float(x[num_start_d50:(num_end_d50 - 1)]))
    if (num_end_particles - num_start_particles == 8):
        particles_file_list = sorted(particles_file_list, key=lambda x: float(x[num_start_particles:num_end_particles]))
    else:
        particles_file_list = sorted(particles_file_list, key=lambda x: float(x[num_start_particles:(num_end_particles - 1)]))
    #print(d50_file_list)
    for x in range(last_timestep_index,len(d50_file_list)):
        rand_num =  random.randint(1, 15)
        print(x)
        d50_nextfile = d50_file_list[x]
        particles_nextfile = particles_file_list[x]
        print(d50_nextfile)
        print(particles_nextfile)
        if (os.path.isfile(d50_nextfile) and os.path.isfile(particles_nextfile)):
            shutil.copy2(d50_nextfile, pasteTo_path)
            shutil.copy2(particles_nextfile, pasteTo_path)
            print("Waiting")
            time.sleep(rand_num)
        else:
            time.sleep(rand_num)
        print("Next file")
    return x


#dummy_PBM(sys.argv[1], np.float(sys.argv[2]), int(sys.argv[3]))
a = dummy_PBM('/home/chai/Documents/git/CyberManufacturing/src/dummy_DEM_PBM/sample_copy',2,1)