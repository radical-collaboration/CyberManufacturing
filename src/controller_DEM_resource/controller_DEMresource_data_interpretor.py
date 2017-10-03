#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 15:32:29 2017

@author: chai
"""

import numpy as np
import controller_DEMresource_data_reader as DEMreader


class controller_DEM_resource_interpretor(object):

    def __init__(self, timestep, types, init_ts):
        self.timestep = timestep
        self.type = types
        self.init_timestep = init_ts
        self.num_of_particles = 0
        self.avg_impacts = 0
        self.tot_part_each_type = np.zeros(self.type)
        self.avg_vel_array = np.zeros(self.type)
        self.collision_matrix = np.zeros((self.type, self.type))

# method to average the velocity of the collision data
    def avg_all_data(self, ts): # here ts is the time step
        obj_data_reader = DEMreader.Controller_DEM_resource_reader(ts, self.type)
        temp_coll_data = np.zeros_like(self.collision_matrix)
        self.num_of_particles = obj_data_reader.number_of_particles
        for x in xrange(0, self.type):
            self.tot_part_each_type[x] = len(obj_data_reader.collision_data_acc_types[str(x+1)])
            temp1 = self.tot_part_each_type[x]
            temp2 = obj_data_reader.collision_data_acc_types[str(x+1)]
            vx = 0.0
            vy = 0.0
            vz = 0.0
            for i in xrange(0, int(temp1)):
                #print(obj_data_reader.collision_data_acc_types[str(x+1)][int(temp1)-1])
                #print(temp2[int(temp1 - 1)]['velocity'][0])
                vz = vz + float(temp2[int(i)]['velocity'][0][2])
                vx = vx + float(temp2[int(i)]['velocity'][0][0])
                vy = vy + float(temp2[int(i)]['velocity'][0][1])
                for j in xrange(0,self.type):
                    temp_coll_data[x][j] += int(temp2[int(i)]['collision'][0][j])
            v_avg = (vx ** 2 + vy ** 2 + vz ** 2) ** (0.5)
            self.avg_vel_array[x] = v_avg / self.tot_part_each_type[x]
        self.collision_matrix = temp_coll_data
        self.avg_impacts = obj_data_reader.number_of_impacts
        print("Recalculated for timestep %d"%ts)

# ------------------------------------------------------------------------------------------------
#
#class controller_DEM_resource_data_init(object):
#    def __init__(self, init_timestep, types):
#        self.init_timestep = init_timestep
#        self.init_impacts = 0
#        self.type = types
#        self.init_collision_matrix = np.zeros((self.type, self.type))
#        self.init_avg_vel = np.zeros(self.type)
#        self.init_num_particles = 0
#
## method to determine the data of the initial time step
#    def init_data_calculations(self):
#        init_data_obj = DEMreader.Controller_DEM_resource_reader(self.init_timestep, self.type)
#        self.init_num_particles = init_data_obj.number_of_particles
#        interpretor_obj = controller_DEM_resource_interpretor(self.init_timestep, self.type)
#        i1 = interpretor_obj.avg_all_data(self.init_timestep)
#        self.init_impacts = init_data_obj.number_of_impacts
#        self.init_avg_vel = interpretor_obj.avg_vel_array
#        self.init_collision_matrix = interpretor_obj.collision_matrix
#        self.init_num_particles = interpretor_obj.num_of_particles
#


# ------------------------------------------------------------------------------------------------

abcd = controller_DEM_resource_interpretor(500000, 16, 200000)
a1 = abcd.avg_all_data(5000000)
#np.set_printoptions(threshold='nan')
print(abcd.collision_matrix)
print(abcd.num_of_particles)
a1 = abcd.avg_all_data(6000000)
print(abcd.collision_matrix)
print(abcd.num_of_particles)