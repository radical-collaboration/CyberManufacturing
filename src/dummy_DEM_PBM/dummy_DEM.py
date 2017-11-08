#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  30 16:19:41 2017

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


def dummy_DEM(pasteTo_path, last_timestep_index, test_number):
    collision_file_list = glob.glob(os.getcwd() + '/test_%d/collision*.atom'%test_number)
    impact_file_list = glob.glob(os.getcwd() + '/test_%d/impact*.atom'%test_number)
    num_start_collision = re.search("collision",collision_file_list[0]).end()
    num_start_impact = re.search("impact",impact_file_list[0]).end()
    num_end_collision = re.search(".atom",collision_file_list[0]).start()
    num_end_impact = re.search(".atom",impact_file_list[0]).start()
    print(impact_file_list)
    if((num_end_collision - num_start_collision) == 7):
        collision_file_list = sorted(collision_file_list, key=lambda x: float(x[(num_start_collision):(num_end_collision + 1)]))
    else:
        collision_file_list = sorted(collision_file_list, key=lambda x: float(x[(num_start_collision):(num_end_collision)]))
    if((num_end_impact - num_start_impact) == 7):
        impact_file_list = sorted(impact_file_list, key=lambda x: float(x[(num_start_impact):(num_end_impact + 1)]))
    else:
        impact_file_list = sorted(impact_file_list, key=lambda x: float(x[(num_start_impact):(num_end_impact)]))
    for x in range(last_timestep_index,len(collision_file_list)):
        rand_num =  random.randint(20, 150)
        print(x)
        collision_nextfile = collision_file_list[x]
        impact_nextfile = impact_file_list[x]
        print(collision_nextfile)
        print(impact_nextfile)
        if (os.path.isfile(collision_nextfile) and os.path.isfile(impact_nextfile)):
            shutil.copy2(collision_nextfile, pasteTo_path)
            shutil.copy2(impact_nextfile, pasteTo_path)
            print("Waiting")
            time.sleep(rand_num)
        else:
            time.sleep(rand_num)
            print("Next file")
            
    return x

#a = dummy_DEM(sys.argv[1], np.float(sys.argv[2]), int(sys.argv[3]))
a = dummy_DEM("/home/chai/Documents/git/CyberManufacturing/src/dummy_DEM_PBM/sample_copy",9000000,1)
