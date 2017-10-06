#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 09:45:37 2017

@author: chai
"""

#this class is for the main execution of the liggghts controller and write data to different files
import numpy as np
import os
import time
import controller_DEMresource_data_reader as DEMreader
import controller_DEMresource_data_interpretor as DEMinter


class controller_DEMresource_main(object):
    def __init__(self, init_timestep, types):
        self.init_timestep = init_timestep
        self.type = types
        self.dump_difference = 50000


    def main(self):
        timestep = self.init_timestep + self.dump_difference
        flag = 0
        dump_avg_vel = open("velocity_average_overtime_starting_%d.txt"%self.init_timestep, "w")
        dump_collisions = open("collisions_overtime_starting_%d.txt"%self.init_timestep, "w")
        dump_impacts = open("impacts_overtime_starting_%d.txt"%self.init_timestep, "w")
        obj_inter = DEMinter.controller_DEM_resource_interpretor(timestep, self.type, self.init_timestep)
        while (flag == 0):
            liggghts_output_files_path = os.getcwd() + '/liggghts_output_files/'
            timestep_collision_file = "collision%d.atom"%timestep
            timestep_impact_file = 'impact%d.atom'%timestep
            collision_file = liggghts_output_files_path + str(timestep_collision_file)
            impact_file = liggghts_output_files_path + str(timestep_impact_file)
            while not(os.path.exists(collision_file) and os.path.exists(impact_file)):
                time.sleep(0.1)
                print("Waiting for file to be printed")
            if (os.path.isfile(collision_file) and os.path.isfile(impact_file)):
                flag = obj_inter.liggghts_data_comparison(timestep)
                a = obj_inter.avg_all_data(timestep)
                dump_avg_vel.write("%d %s \n"%(timestep, str(obj_inter.avg_vel_array)))
                dump_collisions.write("%d %d \n"%(timestep, sum(sum(obj_inter.collision_matrix))))
                dump_impacts.write("%d %d \n"%(timestep, int(obj_inter.avg_impacts)))
                timestep = timestep + self.dump_difference
            else :
                continue
        dump_avg_vel.close()
        dump_collisions.close()
        dump_impacts.close()
        if (flag == 1):
            print("Time to change to PBM and kill DEM")
        elif (flag == 2):
            print("The system is at steady state")


# ------------------------------------------------------------------------------------------------

abcd = controller_DEMresource_main(11500000, 16)
abcd.main()

