#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 09:45:37 2017

@author: Chaitanya Sampat
"""

'''
This class is the main execution file for the DEM controller. It needs initial time and number of type of
particles as inputs. This calls the relevent methods from the interpretor and reader class. It starts 
execution at the inital timestep and then checks for the existence of the file at a predefined timestep
difference (defined as dump_difference). Once the files are found then it calls and the other 2 classes
and checks for the output status from the interpretor class
 ========================================================================================================
 OUTPUT STATUS
 
 0 = continue the execution of the DEM, the % change is less than the defined limit
 1 = kill DEM and execute PBM since the change is  more than 10%
 2 = kill DEM and report that the system has reached a steady state and stop all executions.
 ========================================================================================================
''' 
import numpy as np
import os
import time
import controller_DEMresource_data_reader as DEMreader
import controller_DEMresource_data_interpretor as DEMinter
import json


class controller_DEMresource_main(object):
    def __init__(self, init_timestep, types):
    	# initial timestep and type of partilces are taken as input
        self.init_timestep = init_timestep 
        self.type = types
        self.dump_difference = 50000 # this is the predefined interval after which LIGGGHTS ouputs a dump file. Has to be changed if changed in the LIGGGHTS input file

# ------------------------------------------------------------------------------------------------
    def main(self):
    	# This method executes the necessary commands to run the DEM controller and waits for the command from the interpertor to output certain files that the executioner will act upon.
        timestep = self.init_timestep + self.dump_difference
        flag = 0
        dump_data = {}
        # opening the different files so as to write the calculated average / total data ovetime
        dump_avg_vel = open("velocity_average_overtime_starting_%d.txt"%self.init_timestep, "w") 
        dump_collisions = open("collisions_overtime_starting_%d.txt"%self.init_timestep, "w")
        dump_impacts = open("impacts_overtime_starting_%d.txt"%self.init_timestep, "w")
        obj_inter = DEMinter.controller_DEM_resource_interpretor(timestep, self.type, self.init_timestep)
        # Keeps checking for the existence of the file till one of the criteria for a killing the DEM are not met
        while (flag == 0):
        	# defining the files and path of the files that it needs to search for.
            liggghts_output_files_path = os.getcwd() + '/liggghts_output_files/'
            timestep_collision_file = "collision%d.atom"%timestep
            timestep_impact_file = 'impact%d.atom'%timestep
            collision_file = liggghts_output_files_path + str(timestep_collision_file)
            impact_file = liggghts_output_files_path + str(timestep_impact_file)
            # the execution waits for a seconds everytime it enters the loop till the file do not exist
            while not(os.path.exists(collision_file) and os.path.exists(impact_file)):
                time.sleep(1)
                print("Waiting for file to be printed")
            #once the files are found the it performs the comparison to the initial data files using the interpretor class
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
        	# kill the DEM and start the PBM and also print the input file for the PBM executable
            print("Time to change to PBM and kill DEM")
            tot_part_each_type = np.zeros(16)
            avg_vel_array = np.zeros_like(tot_part_each_type)
            tot_part_each_type = obj_inter.tot_part_each_type
            avg_vel_array = obj_inter.avg_vel_array
            tot_part_each_type.astype(int)
            avg_vel_array.astype(float)
            dump_data['types'] = []
            dump_data['types'].append(self.type)
            dump_data['total number of particles'] = []
            dump_data['total number of particles'].append(int(obj_inter.num_of_particles))
            dump_data['number of particles of each type'] = []
            dump_data['number of particles of each type'].append(list(tot_part_each_type))
            dump_data['average velocity of each type'] = []
            dump_data['average velocity of each type'].append(list(avg_vel_array))
            dump_data['time step'] = []
            dump_data['time step'].append(int(timestep))
            #dump_data['collision_matrix'] = []
            #dump_data['collision_matrix'].append(list(list(obj_inter.collision_matrix)))
            #dump_data.update({'collision_matrix': obj_inter.collision_matrix})
            print(dump_data)
            with open('PBM_input.json' , 'w') as outfile:
            	json.dump(dump_data, outfile, indent=4, separators=(',', ': '))
        elif (flag == 2):
        	# kill the DEM since there has been no change in the number of collisions / impacts / velocity for 5 seconds.
            print("The system is at steady state")


# ------------------------------------------------------------------------------------------------

abcd = controller_DEMresource_main(5000000, 16)
abcd.main()

